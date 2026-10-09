#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12) {
  stop(
    paste(
      "Usage: plot_xci.R",
      "<block_skew.tsv> <summary.tsv> <sensitivity.tsv>",
      "<phase_blocks.tsv> <phase_summary.tsv>",
      "<hp1.bedmethyl.gz> <hp2.bedmethyl.gz>",
      "<cpg_islands.bed> <exclude_beds_csv>",
      "<out_dir> <sample> <bin_bp>"
    )
  )
}

block_file <- args[[1]]
summary_file <- args[[2]]
sensitivity_file <- args[[3]]
phase_file <- args[[4]]
phase_summary_file <- args[[5]]
hp1_file <- args[[6]]
hp2_file <- args[[7]]
cpg_file <- args[[8]]
exclude_arg <- args[[9]]
out_dir <- args[[10]]
sample_id <- args[[11]]
bin_bp <- as.numeric(args[[12]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Thesis-friendly, color-blind-aware palette.
COL_XA <- "#0072B2"
COL_XI <- "#D55E00"
COL_TEAL <- "#009E73"
COL_PURPLE <- "#CC79A7"
COL_GOLD <- "#E69F00"
COL_GREY <- "#9E9E9E"
COL_LIGHT_GREY <- "#D9D9D9"
COL_DARK <- "#333333"

save_plot <- function(plot, name, width, height) {
  prefix <- file.path(out_dir, paste0(sample_id, "_", name))
  ggsave(
    paste0(prefix, ".pdf"),
    plot,
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    paste0(prefix, ".png"),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  ggsave(
    paste0(prefix, ".svg"),
    plot,
    width = width,
    height = height,
    units = "in"
  )
}

theme_thesis <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      axis.title = element_text(size = 11.5),
      axis.text = element_text(size = 10),
      legend.position = "top"
    )
}

first_value <- function(dt, name, default = NA_character_) {
  if (nrow(dt) == 0 || !name %in% names(dt)) {
    return(default)
  }
  as.character(dt[[name]][[1]])
}

blocks <- fread(block_file, na.strings = c("", ".", "NA"))
summary <- fread(summary_file, na.strings = c("", ".", "NA"))
sensitivity <- fread(sensitivity_file, na.strings = c("", ".", "NA"))
phase <- fread(phase_file, na.strings = c("", ".", "NA"))
phase_summary <- fread(
  phase_summary_file,
  na.strings = c("", ".", "NA")
)

global_p <- suppressWarnings(
  as.numeric(first_value(summary, "GLOBAL_FOLDED_SKEW_P", NA))
)
ratio <- first_value(summary, "XCI_MAJOR_MINOR_RATIO", ".")
ci_low <- suppressWarnings(
  as.numeric(
    first_value(
      summary,
      "PROFILE_LIKELIHOOD_CI95_P_LOW",
      NA
    )
  )
)
ci_high <- suppressWarnings(
  as.numeric(
    first_value(
      summary,
      "PROFILE_LIKELIHOOD_CI95_P_HIGH",
      NA
    )
  )
)
lrt_p <- suppressWarnings(
  as.numeric(
    first_value(
      summary,
      "BALANCED_XCI_LRT_BOUNDARY_P",
      NA
    )
  )
)
strong_blocks <- first_value(
  summary,
  "STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_1",
  "."
)

if (nrow(blocks) > 0) {
  numeric_block_cols <- c(
    "FOLDED_BLOCK_SKEW",
    "TRIALS",
    "BLOCK_START",
    "BLOCK_END",
    "H1_XA_PROPORTION_RAW",
    "H1_Xa_SKEW",
    "LOG10_ODDS_H1_XA_VS_H2_XA",
    "PHASE_QC_CONCORDANCE",
    "PHASE_QC_SHARED_SNVS",
    "PHASE_BLOCK_START",
    "PHASE_BLOCK_END"
  )
  for (field in intersect(numeric_block_cols, names(blocks))) {
    blocks[, (field) := as.numeric(get(field))]
  }

  if (!"H1_XA_PROPORTION_RAW" %in% names(blocks) &&
      "H1_Xa_SKEW" %in% names(blocks)) {
    blocks[, H1_XA_PROPORTION_RAW := H1_Xa_SKEW]
  }

  blocks[, BLOCK_START_MB := BLOCK_START / 1e6]
  blocks[, BLOCK_END_MB := BLOCK_END / 1e6]
  blocks[, BLOCK_MID_MB := (BLOCK_START + BLOCK_END) / 2 / 1e6]
}

# -------------------------------------------------------------------------
# 1. Folded block distribution with the global MLE and profile CI.
# -------------------------------------------------------------------------
if (nrow(blocks) > 0) {
  p1 <- ggplot(
    blocks[
      is.finite(FOLDED_BLOCK_SKEW) &
      is.finite(TRIALS) &
      TRIALS > 0
    ]
  ) +
    geom_histogram(
      aes(
        x = FOLDED_BLOCK_SKEW,
        weight = TRIALS
      ),
      binwidth = 0.025,
      boundary = 0,
      fill = COL_TEAL,
      color = "white"
    ) +
    geom_vline(
      xintercept = 0.5,
      linetype = "dotted",
      color = COL_DARK
    )

  if (is.finite(ci_low) && is.finite(ci_high)) {
    p1 <- p1 +
      annotate(
        "rect",
        xmin = ci_low,
        xmax = ci_high,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.14,
        fill = COL_XA
      )
  }

  if (is.finite(global_p)) {
    p1 <- p1 +
      geom_vline(
        xintercept = global_p,
        linetype = "dashed",
        linewidth = 1.0,
        color = COL_XI
      )
  }

  ci_text <- if (is.finite(ci_low) && is.finite(ci_high)) {
    paste0(
      " | 95% profile CI P=",
      sprintf("%.3f", ci_low),
      "–",
      sprintf("%.3f", ci_high)
    )
  } else {
    ""
  }

  lrt_text <- if (is.finite(lrt_p)) {
    paste0(
      " | balanced-XCI LRT p=",
      format(lrt_p, digits = 3, scientific = TRUE)
    )
  } else {
    ""
  }

  p1 <- p1 +
    labs(
      title = paste0(
        sample_id,
        ": X-inactivation skew"
      ),
      subtitle = paste0(
        "Global major:minor estimate ",
        ratio,
        ci_text,
        lrt_text
      ),
      x = "Folded minor-X proportion per block",
      y = "Supporting reads"
    ) +
    coord_cartesian(xlim = c(0, 0.5)) +
    theme_thesis()
} else {
  p1 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No informative X-inactivation blocks"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": X-inactivation skew"
      )
    )
}
save_plot(
  p1,
  "xci_block_skew_distribution",
  10.5,
  6.2
)

# -------------------------------------------------------------------------
# 2. Raw H1-Xa orientation along chrX, showing actual block spans.
# -------------------------------------------------------------------------
if (
  nrow(blocks) > 0 &&
  "H1_XA_PROPORTION_RAW" %in% names(blocks)
) {
  plot_blocks <- blocks[
    is.finite(H1_XA_PROPORTION_RAW) &
    is.finite(BLOCK_START_MB) &
    is.finite(BLOCK_END_MB)
  ]
  plot_blocks[
    ,
    ORIENTATION_DIRECTION := fifelse(
      H1_XA_PROPORTION_RAW >= 0.5,
      "H1 favoured as Xa",
      "H2 favoured as Xa"
    )
  ]

  p2 <- ggplot(
    plot_blocks,
    aes(y = H1_XA_PROPORTION_RAW)
  ) +
    geom_segment(
      aes(
        x = BLOCK_START_MB,
        xend = BLOCK_END_MB,
        yend = H1_XA_PROPORTION_RAW,
        linetype = PHASE_QC_PASS,
        color = ORIENTATION_DIRECTION
      ),
      linewidth = 1.0,
      alpha = 0.7
    ) +
    geom_point(
      aes(
        x = BLOCK_MID_MB,
        size = TRIALS,
        shape = PHASE_QC_PASS,
        color = ORIENTATION_DIRECTION
      ),
      alpha = 0.85
    ) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dashed"
    ) +
    scale_size_continuous(
      name = "Informative reads",
      range = c(2.5, 8)
    ) +
    scale_color_manual(
      name = "Local orientation",
      values = c(
        "H1 favoured as Xa" = COL_XA,
        "H2 favoured as Xa" = COL_XI
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(
        sample_id,
        ": X-inactivation evidence along chrX"
      ),
      subtitle = paste(
        "Horizontal segments show the genomic span of informative CpG",
        "islands in each WhatsHap phase block."
      ),
      x = "chrX position (Mb)",
      y = "H1-as-Xa read proportion",
      linetype = "Phase QC",
      shape = "Phase QC"
    ) +
    theme_thesis()
} else {
  p2 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No informative X-inactivation blocks"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": X-inactivation evidence along chrX"
      )
    )
}
save_plot(
  p2,
  "xci_chrX_block_skew",
  12.5,
  6.0
)

# -------------------------------------------------------------------------
# 3. Signed orientation likelihood along chrX.
# -------------------------------------------------------------------------
if (
  nrow(blocks) > 0 &&
  "LOG10_ODDS_H1_XA_VS_H2_XA" %in% names(blocks)
) {
  odds_blocks <- blocks[
    is.finite(LOG10_ODDS_H1_XA_VS_H2_XA) &
    is.finite(BLOCK_START_MB) &
    is.finite(BLOCK_END_MB)
  ]
  odds_blocks[
    ,
    ORIENTATION_DIRECTION := fifelse(
      LOG10_ODDS_H1_XA_VS_H2_XA >= 0,
      "H1 favoured as Xa",
      "H2 favoured as Xa"
    )
  ]

  p3 <- ggplot(
    odds_blocks,
    aes(y = LOG10_ODDS_H1_XA_VS_H2_XA)
  ) +
    geom_segment(
      aes(
        x = BLOCK_START_MB,
        xend = BLOCK_END_MB,
        yend = LOG10_ODDS_H1_XA_VS_H2_XA,
        linetype = PHASE_QC_PASS,
        color = ORIENTATION_DIRECTION
      ),
      linewidth = 1.0,
      alpha = 0.7
    ) +
    geom_point(
      aes(
        x = BLOCK_MID_MB,
        size = TRIALS,
        shape = PHASE_QC_PASS,
        color = ORIENTATION_DIRECTION
      ),
      alpha = 0.85
    ) +
    geom_hline(yintercept = 0) +
    geom_hline(
      yintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = c(-2, 2),
      linetype = "dotted"
    ) +
    scale_size_continuous(
      name = "Informative reads",
      range = c(2.5, 8)
    ) +
    scale_color_manual(
      name = "Local orientation",
      values = c(
        "H1 favoured as Xa" = COL_XA,
        "H2 favoured as Xa" = COL_XI
      )
    ) +
    labs(
      title = paste0(
        sample_id,
        ": confidence in X-inactivation orientation"
      ),
      subtitle = paste0(
        "|log10 odds| ≥1 = at least 10:1 support; ≥2 = at least 100:1. ",
        "Strong blocks: ",
        strong_blocks
      ),
      x = "chrX position (Mb)",
      y = "log10 LR: H1 Xa vs H2 Xa",
      linetype = "Phase QC",
      shape = "Phase QC"
    ) +
    theme_thesis()
} else {
  p3 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No block-level orientation odds available"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": confidence in X-inactivation orientation"
      )
    )
}
save_plot(
  p3,
  "xci_orientation_log_odds",
  12.5,
  6.0
)

# -------------------------------------------------------------------------
# Helpers for CpG-island-only, block-oriented Xa/Xi methylation.
# -------------------------------------------------------------------------
read_cpg_bed <- function(path) {
  dt <- fread(
    path,
    header = FALSE,
    fill = TRUE,
    showProgress = FALSE
  )
  dt <- dt[
    !grepl("^#", as.character(V1))
  ]
  if (ncol(dt) < 3) {
    stop("CpG-island BED needs at least 3 columns")
  }
  dt <- dt[, .(
    chrom = as.character(V1),
    cpg_start = as.numeric(V2),
    cpg_end = as.numeric(V3)
  )]
  dt <- dt[
    chrom == "chrX" &
    is.finite(cpg_start) &
    is.finite(cpg_end) &
    cpg_end > cpg_start
  ]
  unique(dt)
}

read_exclusions <- function(csv) {
  if (is.na(csv) || csv == "" || csv == "NONE") {
    return(
      data.table(
        chrom = character(),
        ex_start = numeric(),
        ex_end = numeric()
      )
    )
  }

  paths <- strsplit(csv, ",", fixed = TRUE)[[1]]
  out <- list()

  for (path in paths[nzchar(paths)]) {
    if (!file.exists(path)) {
      stop(
        paste(
          "XCI exclusion BED does not exist:",
          path
        )
      )
    }

    dt <- fread(
      path,
      header = FALSE,
      fill = TRUE,
      showProgress = FALSE
    )
    dt <- dt[
      !grepl("^#", as.character(V1))
    ]
    if (ncol(dt) < 3) {
      stop(
        paste(
          "XCI exclusion BED needs at least 3 columns:",
          path
        )
      )
    }

    out[[length(out) + 1]] <- dt[, .(
      chrom = as.character(V1),
      ex_start = as.numeric(V2),
      ex_end = as.numeric(V3)
    )]
  }

  rbindlist(out, fill = TRUE)[
    chrom == "chrX" &
    is.finite(ex_start) &
    is.finite(ex_end) &
    ex_end > ex_start
  ]
}

mask_cpgs <- function(cpgs, exclusions) {
  if (nrow(cpgs) == 0 || nrow(exclusions) == 0) {
    return(cpgs)
  }

  x <- copy(cpgs)
  x[, row_id := .I]

  setkey(exclusions, chrom, ex_start, ex_end)
  hit <- foverlaps(
    x,
    exclusions,
    by.x = c("chrom", "cpg_start", "cpg_end"),
    by.y = c("chrom", "ex_start", "ex_end"),
    type = "any",
    nomatch = 0L
  )

  remove_ids <- unique(hit$row_id)
  x <- x[!row_id %in% remove_ids]
  x[, row_id := NULL]
  x
}

read_bedmethyl <- function(path, haplotype) {
  cmd <- paste("gzip -dc", shQuote(path))
  dt <- fread(
    cmd = cmd,
    header = FALSE,
    select = 1:11,
    showProgress = FALSE
  )
  if (ncol(dt) < 11) {
    stop(
      paste(
        "bedMethyl has fewer than 11 columns:",
        path
      )
    )
  }

  setnames(
    dt,
    c(
      "chrom", "start", "end", "name", "score", "strand",
      "thick_start", "thick_end", "color", "coverage",
      "methylation"
    )
  )

  dt <- dt[
    chrom == "chrX"
  ]
  dt[, start := as.numeric(start)]
  dt[, end := as.numeric(end)]
  dt[, coverage := as.numeric(coverage)]
  dt[, methylation := as.numeric(methylation)]

  dt <- dt[
    is.finite(start) &
    is.finite(end) &
    is.finite(coverage) &
    is.finite(methylation) &
    coverage > 0
  ]

  dt[, haplotype := haplotype]
  dt[, row_id := paste0(haplotype, "_", .I)]
  dt
}

cpgs <- mask_cpgs(
  read_cpg_bed(cpg_file),
  read_exclusions(exclude_arg)
)

orientable <- copy(blocks)
if (nrow(orientable) > 0) {
  orientable <- orientable[
    ORIENTATION_USABLE == "YES" &
    PREFERRED_XA_HAPLOTYPE %in% c("H1", "H2") &
    is.finite(BLOCK_START) &
    is.finite(BLOCK_END)
  ]
  orientable[, ABS_LOG10_ODDS := abs(
    LOG10_ODDS_H1_XA_VS_H2_XA
  )]
  orientable[, chrom := "chrX"]
  orientable[
    ,
    block_start := fifelse(
      is.finite(PHASE_BLOCK_START),
      PHASE_BLOCK_START,
      BLOCK_START
    )
  ]
  orientable[
    ,
    block_end := fifelse(
      is.finite(PHASE_BLOCK_END),
      PHASE_BLOCK_END,
      BLOCK_END
    )
  ]
}

hp1 <- read_bedmethyl(hp1_file, "HP1")
hp2 <- read_bedmethyl(hp2_file, "HP2")
meth_raw <- rbindlist(
  list(hp1, hp2),
  fill = TRUE
)

# CpG-island restriction.
if (nrow(meth_raw) > 0 && nrow(cpgs) > 0) {
  setkey(cpgs, chrom, cpg_start, cpg_end)
  meth_cgi <- foverlaps(
    meth_raw,
    cpgs,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "cpg_start", "cpg_end"),
    type = "any",
    nomatch = 0L
  )
  meth_cgi <- unique(
    meth_cgi,
    by = "row_id"
  )
} else {
  meth_cgi <- data.table()
}

# Assign each methylation record to a confidently oriented WhatsHap block.
oriented_meth <- data.table()
if (
  nrow(meth_cgi) > 0 &&
  nrow(orientable) > 0
) {
  setkey(
    orientable,
    chrom,
    block_start,
    block_end
  )

  oriented_meth <- foverlaps(
    meth_cgi,
    orientable,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "block_start", "block_end"),
    type = "within",
    nomatch = 0L
  )

  # If a site is covered by more than one candidate block, retain the
  # orientation with the strongest absolute likelihood ratio.
  setorder(
    oriented_meth,
    row_id,
    -ABS_LOG10_ODDS
  )
  oriented_meth <- unique(
    oriented_meth,
    by = "row_id"
  )

  oriented_meth[, X_STATE := fifelse(
    (
      haplotype == "HP1" &
      PREFERRED_XA_HAPLOTYPE == "H1"
    ) |
    (
      haplotype == "HP2" &
      PREFERRED_XA_HAPLOTYPE == "H2"
    ),
    "Xa",
    "Xi"
  )]

  oriented_meth[, bin_start := floor(start / bin_bp) * bin_bp]
  oriented_meth[
    ,
    weighted_methylation := methylation * coverage
  ]
}

if (nrow(oriented_meth) > 0) {
  meth_summary <- oriented_meth[
    ,
    .(
      coverage_sum = sum(coverage, na.rm = TRUE),
      weighted_sum = sum(
        weighted_methylation,
        na.rm = TRUE
      ),
      CpG_records = .N,
      phase_blocks = uniqueN(PS)
    ),
    by = .(
      bin_start,
      X_STATE
    )
  ]
  meth_summary[
    ,
    methylation_percent := weighted_sum / coverage_sum
  ]
  meth_summary[
    ,
    position_mb := (bin_start + bin_bp / 2) / 1e6
  ]

  fwrite(
    meth_summary,
    file.path(
      out_dir,
      paste0(
        sample_id,
        "_xci_oriented_CGI_methylation.tsv"
      )
    ),
    sep = "\t",
    na = "."
  )

  p4 <- ggplot(
    meth_summary[
      is.finite(methylation_percent)
    ],
    aes(
      position_mb,
      methylation_percent,
      group = X_STATE,
      linetype = X_STATE,
      color = X_STATE
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_color_manual(
      name = "Oriented X state",
      values = c(
        "Xa" = COL_XA,
        "Xi" = COL_XI
      )
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      title = paste0(
        sample_id,
        ": Xa/Xi CpG-island methylation"
      ),
      subtitle = paste(
        "HP1/HP2 were re-oriented within phase-QC-passed blocks",
        "with |log10 odds| ≥ 1; only non-masked chrX CpG islands are shown."
      ),
      x = "chrX position (Mb)",
      y = "CpG-island 5mC (%)",
      linetype = "Oriented X state",
      color = "Oriented X state"
    ) +
    theme_thesis()

  wide <- dcast(
    meth_summary,
    bin_start + position_mb ~ X_STATE,
    value.var = "methylation_percent"
  )

  if (
    "Xa" %in% names(wide) &&
    "Xi" %in% names(wide)
  ) {
    wide[, XI_MINUS_XA := Xi - Xa]
    wide[
      ,
      DELTA_DIRECTION := fifelse(
        XI_MINUS_XA >= 0,
        "Xi > Xa",
        "Xa > Xi"
      )
    ]

    fwrite(
      wide,
      file.path(
        out_dir,
        paste0(
          sample_id,
          "_xci_methylation_delta.tsv"
        )
      ),
      sep = "\t",
      na = "."
    )

    p4_delta <- ggplot(
      wide[
        is.finite(XI_MINUS_XA)
      ],
      aes(
        position_mb,
        XI_MINUS_XA,
        color = DELTA_DIRECTION
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.2) +
      scale_color_manual(
        name = "Methylation difference",
        values = c(
          "Xi > Xa" = COL_XI,
          "Xa > Xi" = COL_XA
        )
      ) +
      labs(
        title = paste0(
          sample_id,
          ": Xi − Xa CpG-island methylation"
        ),
        subtitle = paste(
          "Positive values indicate higher CpG-island 5mC on the",
          "locally oriented inactive X."
        ),
        x = "chrX position (Mb)",
        y = "Xi − Xa 5mC (percentage points)"
      ) +
      theme_thesis()
  } else {
    p4_delta <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = "No bins contain both confidently oriented Xa and Xi methylation"
      ) +
      theme_void() +
      labs(
        title = paste0(
          sample_id,
          ": Xi − Xa CpG-island methylation"
        )
      )
  }
} else {
  p4 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste(
        "No CpG methylation records remained after",
        "block orientation and phase-QC filtering"
      )
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": Xa/Xi CpG-island methylation"
      )
    )

  p4_delta <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No oriented CpG-island methylation available"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": Xi − Xa CpG-island methylation"
      )
    )
}

save_plot(
  p4,
  "xci_haplotype_methylation",
  12.5,
  6.0
)
save_plot(
  p4_delta,
  "xci_methylation_delta",
  12.5,
  5.5
)

# -------------------------------------------------------------------------
# 5. WhatsHap vs LongPhase concordance: show all overlap pairs and n.
# -------------------------------------------------------------------------
if (nrow(phase) > 0) {
  phase[, ORIENTATION_CONCORDANCE := as.numeric(
    ORIENTATION_CONCORDANCE
  )]
  phase[, SHARED_PHASED_SNVS := as.numeric(
    SHARED_PHASED_SNVS
  )]

  phase <- phase[
    order(
      SHARED_PHASED_SNVS,
      ORIENTATION_CONCORDANCE
    )
  ]
  phase[
    ,
    block := paste0(
      "W",
      WHATSHAP_PS,
      "/L",
      LONGPHASE_PS
    )
  ]
  phase[
    ,
    block := factor(
      block,
      levels = block
    )
  ]

  overall <- first_value(
    phase_summary,
    "FLIP_TOLERANT_PHASE_CONCORDANCE",
    "."
  )
  phase[
    ,
    CONCORDANCE_CLASS := fifelse(
      ORIENTATION_CONCORDANCE >= 0.90,
      "≥ 0.90 concordance",
      "< 0.90 concordance"
    )
  ]

  p5 <- ggplot(
    phase,
    aes(
      block,
      ORIENTATION_CONCORDANCE,
      fill = CONCORDANCE_CLASS
    )
  ) +
    geom_col() +
    scale_fill_manual(
      name = "Phase QC",
      values = c(
        "≥ 0.90 concordance" = COL_TEAL,
        "< 0.90 concordance" = COL_GREY
      )
    ) +
    geom_text(
      aes(
        label = paste0(
          "n=",
          SHARED_PHASED_SNVS
        )
      ),
      hjust = -0.08,
      size = 3.0
    ) +
    coord_flip(
      ylim = c(0.5, 1.07),
      clip = "off"
    ) +
    labs(
      title = paste0(
        sample_id,
        ": WhatsHap and LongPhase agreement"
      ),
      subtitle = paste0(
        "Overall flip-tolerant shared-SNV concordance: ",
        overall,
        ". All overlapping block pairs are shown."
      ),
      x = "Overlapping phase-block pair",
      y = "Flip-tolerant phase concordance"
    ) +
    theme_thesis()

  phase_height <- max(
    6.5,
    2.5 + 0.22 * nrow(phase)
  )
} else {
  p5 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No shared chrX phase blocks"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": WhatsHap and LongPhase agreement"
      )
    )
  phase_height <- 6.5
}
save_plot(
  p5,
  "xci_phase_concordance",
  11.5,
  phase_height
)

# -------------------------------------------------------------------------
# 6. Sensitivity to read threshold and phase-QC filtering.
# -------------------------------------------------------------------------
if (nrow(sensitivity) > 0) {
  for (field in c(
    "MIN_READS_PER_BLOCK",
    "ALL_BLOCKS_FOLDED_MINOR_X_P",
    "PHASE_QC_FOLDED_MINOR_X_P"
  )) {
    if (field %in% names(sensitivity)) {
      sensitivity[
        ,
        (field) := as.numeric(get(field))
      ]
    }
  }

  sens_long <- melt(
    sensitivity,
    id.vars = "MIN_READS_PER_BLOCK",
    measure.vars = c(
      "ALL_BLOCKS_FOLDED_MINOR_X_P",
      "PHASE_QC_FOLDED_MINOR_X_P"
    ),
    variable.name = "analysis",
    value.name = "minor_p"
  )
  sens_long[
    ,
    analysis := fifelse(
      analysis == "ALL_BLOCKS_FOLDED_MINOR_X_P",
      "All read-qualified blocks",
      "Phase-QC-passed blocks"
    )
  ]

  p6 <- ggplot(
    sens_long[
      is.finite(minor_p)
    ],
    aes(
      MIN_READS_PER_BLOCK,
      minor_p,
      group = analysis,
      linetype = analysis,
      shape = analysis,
      color = analysis
    )
  ) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dotted"
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8) +
    scale_color_manual(
      name = "Block set",
      values = c(
        "All read-qualified blocks" = COL_XA,
        "Phase-QC-passed blocks" = COL_TEAL
      )
    ) +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(
      title = paste0(
        sample_id,
        ": XCI estimate sensitivity"
      ),
      subtitle = paste(
        "Primary result remains the folded-binomial MLE using the configured",
        "minimum-read threshold; phase filtering is shown only as sensitivity."
      ),
      x = "Minimum informative reads per block",
      y = "Folded minor-X P",
      linetype = "Block set",
      shape = "Block set",
      color = "Block set"
    ) +
    theme_thesis()
} else {
  p6 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No XCI sensitivity results available"
    ) +
    theme_void() +
    labs(
      title = paste0(
        sample_id,
        ": XCI estimate sensitivity"
      )
    )
}
save_plot(
  p6,
  "xci_sensitivity",
  10.5,
  5.8
)

manifest <- data.table(
  sample = sample_id,
  XCI_MAJOR_MINOR_RATIO = ratio,
  GLOBAL_FOLDED_SKEW_P = ifelse(
    is.finite(global_p),
    global_p,
    NA_real_
  ),
  PROFILE_LIKELIHOOD_CI95_P_LOW = ci_low,
  PROFILE_LIKELIHOOD_CI95_P_HIGH = ci_high,
  BALANCED_XCI_LRT_BOUNDARY_P = lrt_p,
  plot_directory = out_dir
)

fwrite(
  manifest,
  file.path(
    out_dir,
    paste0(
      sample_id,
      "_xci_plot_manifest.tsv"
    )
  ),
  sep = "\t",
  na = "."
)

message(
  "[OK] X-inactivation R plots written to ",
  out_dir
)
