#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(NanoMethViz)
  library(tidyverse)
  library(doParallel)
  library(foreach)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop(
    paste(
      "Usage: xci_cluster_methylation.R",
      "<sample> <bam> <cpg_bed> <haplotags.tsv.gz>",
      "<clustered.tsv.gz> <block_skew.tsv.gz>",
      "<min_cluster_reads> <threads> <xist_promoter_bed_or_NONE>",
      "<comma_separated_exclude_beds_or_NONE>"
    )
  )
}

sample_id <- args[[1]]
bam_path <- args[[2]]
cpg_bed_path <- args[[3]]
haplotag_path <- args[[4]]
clustered_out <- args[[5]]
block_out <- args[[6]]
min_pts <- as.integer(args[[7]])
nthreads <- as.integer(args[[8]])
xist_bed_path <- args[[9]]
exclude_bed_arg <- args[[10]]

cpgs <- read_tsv(
  cpg_bed_path,
  col_names = FALSE,
  comment = "#",
  show_col_types = FALSE
)

if (ncol(cpgs) < 3) {
  stop("CpG-island BED requires at least chromosome, start and end columns")
}

cpgs <- cpgs %>%
  transmute(
    chr = as.character(X1),
    start = as.integer(X2),
    end = as.integer(X3)
  ) %>%
  filter(chr == "chrX", end > start) %>%
  distinct()

if (nrow(cpgs) == 0) {
  stop("No chrX intervals were found in the configured CpG-island BED")
}

read_exclude_beds <- function(arg) {
  if (is.na(arg) || arg == "" || arg == "NONE") {
    return(tibble(chr = character(), start = integer(), end = integer(), label = character()))
  }

  paths <- strsplit(arg, ",", fixed = TRUE)[[1]]
  paths <- paths[nzchar(paths)]
  out <- list()

  for (path in paths) {
    if (!file.exists(path)) {
      stop(paste("Configured XCI exclusion BED does not exist:", path))
    }

    raw <- read_tsv(
      path,
      col_names = FALSE,
      comment = "#",
      show_col_types = FALSE
    )
    if (ncol(raw) < 3) {
      stop(paste("XCI exclusion BED needs at least 3 columns:", path))
    }

    label <- if (ncol(raw) >= 4) as.character(raw[[4]]) else basename(path)
    out[[length(out) + 1]] <- tibble(
      chr = as.character(raw[[1]]),
      start = as.integer(raw[[2]]),
      end = as.integer(raw[[3]]),
      label = label
    )
  }

  bind_rows(out) %>%
    filter(chr == "chrX", end > start)
}

exclude_bed <- read_exclude_beds(exclude_bed_arg)

overlaps_excluded <- function(chr, start, end) {
  if (nrow(exclude_bed) == 0) {
    return(FALSE)
  }
  any(
    exclude_bed$chr == chr &
      exclude_bed$end > start &
      exclude_bed$start < end
  )
}

before_exclusion <- nrow(cpgs)
cpgs <- cpgs %>%
  rowwise() %>%
  filter(!overlaps_excluded(chr, start, end)) %>%
  ungroup()

message(
  sprintf(
    "[XCI] CpG islands retained after PAR/XIST/escape masking: %d/%d",
    nrow(cpgs),
    before_exclusion
  )
)

if (nrow(cpgs) == 0) {
  stop("All chrX CpG islands were removed by the configured XCI exclusion masks")
}

xist_bed <- tibble(chr = character(), start = integer(), end = integer())
if (
  !is.na(xist_bed_path) &&
  xist_bed_path != "NONE" &&
  file.exists(xist_bed_path)
) {
  xraw <- read_tsv(
    xist_bed_path,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE
  )
  if (ncol(xraw) >= 3) {
    xist_bed <- xraw %>%
      transmute(
        chr = as.character(X1),
        start = as.integer(X2),
        end = as.integer(X3)
      ) %>%
      filter(chr == "chrX", end > start)
  }
}

is_xist_interval <- function(chr, start, end) {
  if (nrow(xist_bed) == 0) {
    return(FALSE)
  }
  any(
    xist_bed$chr == chr &
      xist_bed$end > start &
      xist_bed$start < end
  )
}

mbr <- ModBamResult(
  methy = ModBamFiles(
    samples = sample_id,
    paths = bam_path
  ),
  samples = data.frame(
    sample = sample_id,
    group = 1
  ),
  # Dorado BAMs may contain both 5mC (MM code m) and 5hmC (MM code h).
  # XCI clustering follows the SkewX method and uses CpG 5mC specifically.
  # Setting this explicitly avoids ambiguity in multi-modification modBAMs.
  mod_code = "m"
)

cluster_one <- function(row_index) {
  region <- cpgs[row_index, ]

  clustered <- tryCatch(
    NanoMethViz:::cluster_reads(
      mbr,
      region$chr,
      region$start,
      region$end,
      min_pts = min_pts
    ),
    error = function(e) {
      message(
        sprintf(
          "Skipping %s:%s-%s: %s",
          region$chr,
          region$start,
          region$end,
          e$message
        )
      )
      NULL
    }
  )

  if (is.null(clustered) || nrow(clustered) == 0) {
    return(NULL)
  }

  required <- c("read_name", "cluster_id", "mean")
  missing <- setdiff(required, colnames(clustered))
  if (length(missing) > 0) {
    stop(
      paste(
        "NanoMethViz cluster_reads output is missing columns:",
        paste(missing, collapse = ", ")
      )
    )
  }

  clustered <- clustered %>%
    mutate(
      CGI_id = paste0(
        region$chr, ":", region$start, "-", region$end
      ),
      chr = region$chr,
      start = region$start,
      end = region$end,
      cluster_id = as.character(cluster_id)
    )

  valid_cluster_ids <- sort(
    unique(
      clustered$cluster_id[
        !is.na(clustered$cluster_id) &
          clustered$cluster_id != ""
      ]
    )
  )

  clustered <- clustered %>%
    group_by(cluster_id) %>%
    mutate(
      avg_cluster_methylation = mean(mean, na.rm = TRUE)
    ) %>%
    ungroup()

  xist_reverse <- is_xist_interval(
    region$chr,
    region$start,
    region$end
  )

  if (length(valid_cluster_ids) == 2) {
    cluster_means <- clustered %>%
      filter(cluster_id %in% valid_cluster_ids) %>%
      group_by(cluster_id) %>%
      summarise(
        avg = mean(avg_cluster_methylation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(avg)

    low_cluster <- cluster_means$cluster_id[[1]]
    high_cluster <- cluster_means$cluster_id[[2]]

    clustered <- clustered %>%
      mutate(
        assigned_X = case_when(
          !xist_reverse & cluster_id == low_cluster ~ "Xa",
          !xist_reverse & cluster_id == high_cluster ~ "Xi",
          xist_reverse & cluster_id == low_cluster ~ "Xi",
          xist_reverse & cluster_id == high_cluster ~ "Xa",
          TRUE ~ NA_character_
        ),
        XCI_assignment_rule = ifelse(
          xist_reverse,
          "XIST_PROMOTER_REVERSED",
          "STANDARD_LOW_Xa_HIGH_Xi"
        )
      )
  } else {
    clustered <- clustered %>%
      mutate(
        assigned_X = NA_character_,
        XCI_assignment_rule = "NOT_EXACTLY_TWO_CLUSTERS"
      )
  }

  clustered
}

cl <- makeCluster(max(1, nthreads))
registerDoParallel(cl)

clustered_list <- foreach(
  idx = seq_len(nrow(cpgs)),
  .packages = c("NanoMethViz", "tidyverse")
) %dopar% {
  cluster_one(idx)
}

stopCluster(cl)

clustered_reads <- bind_rows(clustered_list)

cluster_fields <- c(
  "read_name", "cluster_id", "mean",
  "avg_cluster_methylation", "assigned_X",
  "XCI_assignment_rule", "CGI_id",
  "chr", "start", "end"
)

if (nrow(clustered_reads) == 0) {
  clustered_reads <- as_tibble(
    setNames(
      replicate(
        length(cluster_fields),
        character(),
        simplify = FALSE
      ),
      cluster_fields
    )
  )
}

write_tsv(clustered_reads, clustered_out)

haplotags <- read_tsv(
  haplotag_path,
  show_col_types = FALSE,
  col_types = cols(
    read_name = col_character(),
    HP = col_integer(),
    PS = col_character(),
    .default = col_guess()
  )
)

informative <- clustered_reads %>%
  filter(assigned_X %in% c("Xa", "Xi")) %>%
  inner_join(
    haplotags %>% select(read_name, HP, PS),
    by = "read_name"
  ) %>%
  filter(
    HP %in% c(1L, 2L),
    !is.na(PS),
    PS != "."
  )

if (nrow(informative) == 0) {
  empty_block <- tibble(
    PS = character(),
    BLOCK_START = integer(),
    BLOCK_END = integer(),
    INFORMATIVE_CPG_ISLANDS = integer(),
    UNIQUE_READS = integer(),
    MULTI_CGI_READS = integer(),
    DISCORDANT_MULTI_CGI_READS = integer(),
    CONSENSUS_READS = integer(),
    H1_Xa = integer(),
    H1_Xi = integer(),
    H2_Xa = integer(),
    H2_Xi = integer(),
    H1_Xa_SKEW = double()
  )
  write_tsv(empty_block, block_out)
  quit(save = "no", status = 0)
}

block_meta <- informative %>%
  group_by(PS) %>%
  summarise(
    BLOCK_START = min(start, na.rm = TRUE),
    BLOCK_END = max(end, na.rm = TRUE),
    INFORMATIVE_CPG_ISLANDS = n_distinct(CGI_id),
    UNIQUE_READS = n_distinct(read_name),
    .groups = "drop"
  )

# Match the published SkewX rule that each read contributes once per
# haplotype block, but do not resolve multi-island conflicts by string order.
# If the same read is assigned Xa at one informative CGI and Xi at another,
# its Xa/Xi state is discordant and the read is excluded from the skew count.
per_read_consensus <- informative %>%
  group_by(read_name, HP, PS) %>%
  summarise(
    N_CGI_ASSIGNMENTS = n_distinct(CGI_id),
    N_X_STATES = n_distinct(assigned_X),
    CONSENSUS_X = ifelse(
      N_X_STATES == 1,
      first(assigned_X),
      NA_character_
    ),
    .groups = "drop"
  )

read_qc <- per_read_consensus %>%
  group_by(PS) %>%
  summarise(
    MULTI_CGI_READS = sum(N_CGI_ASSIGNMENTS > 1),
    DISCORDANT_MULTI_CGI_READS = sum(
      N_CGI_ASSIGNMENTS > 1 & is.na(CONSENSUS_X)
    ),
    CONSENSUS_READS = sum(!is.na(CONSENSUS_X)),
    .groups = "drop"
  )

informative_once <- per_read_consensus %>%
  filter(!is.na(CONSENSUS_X)) %>%
  transmute(
    read_name = read_name,
    HP = HP,
    PS = PS,
    assigned_X = CONSENSUS_X
  )

counts <- informative_once %>%
  count(PS, assigned_X, HP, name = "counts") %>%
  mutate(
    combination = case_when(
      assigned_X == "Xa" & HP == 1 ~ "H1_Xa",
      assigned_X == "Xi" & HP == 1 ~ "H1_Xi",
      assigned_X == "Xa" & HP == 2 ~ "H2_Xa",
      assigned_X == "Xi" & HP == 2 ~ "H2_Xi",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(combination)) %>%
  select(PS, combination, counts) %>%
  pivot_wider(
    names_from = combination,
    values_from = counts,
    values_fill = 0
  )

for (field in c("H1_Xa", "H1_Xi", "H2_Xa", "H2_Xi")) {
  if (!field %in% colnames(counts)) {
    counts[[field]] <- 0L
  }
}

block_skew <- counts %>%
  left_join(block_meta, by = "PS") %>%
  left_join(read_qc, by = "PS") %>%
  mutate(
    TRIALS = H1_Xa + H1_Xi + H2_Xa + H2_Xi,
    H1_Xa_SKEW = ifelse(
      TRIALS > 0,
      (H1_Xa + H2_Xi) / TRIALS,
      NA_real_
    )
  ) %>%
  select(
    PS,
    BLOCK_START,
    BLOCK_END,
    INFORMATIVE_CPG_ISLANDS,
    UNIQUE_READS,
    MULTI_CGI_READS,
    DISCORDANT_MULTI_CGI_READS,
    CONSENSUS_READS,
    H1_Xa,
    H1_Xi,
    H2_Xa,
    H2_Xi,
    H1_Xa_SKEW
  ) %>%
  arrange(BLOCK_START, BLOCK_END)

write_tsv(block_skew, block_out)

message(
  sprintf(
    "[XCI] per-read consensus: retained=%d discordant_multi_CGI=%d",
    sum(block_skew$CONSENSUS_READS, na.rm = TRUE),
    sum(block_skew$DISCORDANT_MULTI_CGI_READS, na.rm = TRUE)
  )
)
