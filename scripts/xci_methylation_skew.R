#!/usr/bin/env Rscript

# Allele-specific X-chromosome inactivation analysis from native ONT reads.
#
# Methodological basis:
# - cluster read-level CpG-island methylation into two epiallele groups;
# - assign lower methylation to Xa and higher methylation to Xi;
# - reverse that interpretation for configured XIST-promoter CpG islands;
# - combine Xa/Xi labels with WhatsHap HP/PS read tags;
# - estimate block-wise and sample-wide XCI skew with a folded-binomial MLE.
#
# This is a research measurement of XCI skew. It is not a pathogenicity score.

suppressPackageStartupMessages({
    library(NanoMethViz)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(foreach)
    library(doParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
    stop(
        paste(
            "Usage: xci_methylation_skew.R SAMPLE BAM CPG_BED HAPLOTAGS",
            "CLUSTERED_OUT BLOCK_OUT SUMMARY_OUT MIN_CLUSTER_READS",
            "MIN_BLOCK_READS THREADS [XIST_PROMOTER_BED_OR_DOT]"
        )
    )
}

sample_id <- args[[1]]
bam <- args[[2]]
cpg_bed_path <- args[[3]]
haplotag_path <- args[[4]]
clustered_out <- args[[5]]
block_out <- args[[6]]
summary_out <- args[[7]]
min_cluster_reads <- as.integer(args[[8]])
min_block_reads <- as.integer(args[[9]])
threads <- as.integer(args[[10]])
xist_bed_path <- ifelse(length(args) >= 11, args[[11]], ".")

read_bed3 <- function(path) {
    x <- suppressMessages(
        read_tsv(
            path,
            col_names = FALSE,
            comment = "#",
            show_col_types = FALSE,
            progress = FALSE
        )
    )
    if (ncol(x) < 3) {
        stop(paste("BED needs at least 3 columns:", path))
    }
    x <- x[, 1:3]
    names(x) <- c("CHROM", "START", "END")
    x %>%
        mutate(
            CHROM = as.character(CHROM),
            START = as.integer(START),
            END = as.integer(END)
        ) %>%
        filter(!is.na(START), !is.na(END), END > START)
}

cpg_bed <- read_bed3(cpg_bed_path)
if (nrow(cpg_bed) == 0) {
    stop("CpG-island BED contains no usable intervals.")
}

xist_bed <- NULL
if (!is.na(xist_bed_path) && xist_bed_path != "." && file.exists(xist_bed_path)) {
    xist_bed <- read_bed3(xist_bed_path)
}

haplotags <- suppressMessages(
    read_tsv(
        haplotag_path,
        show_col_types = FALSE,
        progress = FALSE
    )
)
required_hap <- c("read_name", "HP", "PS")
if (!all(required_hap %in% names(haplotags))) {
    stop(
        paste(
            "Haplotype table is missing:",
            paste(setdiff(required_hap, names(haplotags)), collapse = ", ")
        )
    )
}
haplotags <- haplotags %>%
    transmute(
        read_name = as.character(read_name),
        HP = as.integer(HP),
        PS = as.character(PS)
    ) %>%
    filter(HP %in% c(1L, 2L), !is.na(PS), PS != ".")

mbr <- ModBamResult(
    methy = ModBamFiles(samples = sample_id, paths = bam),
    samples = data.frame(sample = sample_id, group = 1)
)

overlaps_xist <- function(chrom, start, end) {
    if (is.null(xist_bed) || nrow(xist_bed) == 0) {
        return(FALSE)
    }
    any(
        xist_bed$CHROM == chrom &
        xist_bed$START < end &
        xist_bed$END > start
    )
}

process_cpg <- function(i) {
    current <- cpg_bed[i, ]

    clustered <- tryCatch(
        NanoMethViz:::cluster_reads(
            mbr,
            current$CHROM,
            current$START,
            current$END,
            min_pts = min_cluster_reads
        ),
        error = function(e) NULL
    )

    if (is.null(clustered) || nrow(clustered) == 0) {
        return(NULL)
    }
    if (!all(c("read_name", "cluster_id", "mean") %in% names(clustered))) {
        return(NULL)
    }

    clustered <- clustered %>%
        mutate(
            cluster_id = as.character(cluster_id),
            CGI_ID = paste0(
                current$CHROM, ":", current$START, "-", current$END
            ),
            CHROM = current$CHROM,
            START = current$START,
            END = current$END
        )

    valid_clusters <- sort(
        unique(clustered$cluster_id[!is.na(clustered$cluster_id)])
    )

    clustered$assigned_X <- NA_character_
    clustered$XIST_LABEL_REVERSED <- "NO"

    if (length(valid_clusters) == 2) {
        cluster_means <- clustered %>%
            filter(cluster_id %in% valid_clusters) %>%
            group_by(cluster_id) %>%
            summarise(
                CLUSTER_MEAN_METHYLATION = mean(mean, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            arrange(CLUSTER_MEAN_METHYLATION)

        low_cluster <- cluster_means$cluster_id[[1]]
        high_cluster <- cluster_means$cluster_id[[2]]
        is_xist <- overlaps_xist(
            current$CHROM,
            current$START,
            current$END
        )

        if (is_xist) {
            # XIST promoter methylation has the opposite Xa/Xi relationship.
            clustered$assigned_X <- case_when(
                clustered$cluster_id == low_cluster ~ "Xi",
                clustered$cluster_id == high_cluster ~ "Xa",
                TRUE ~ NA_character_
            )
            clustered$XIST_LABEL_REVERSED <- "YES"
        } else {
            clustered$assigned_X <- case_when(
                clustered$cluster_id == low_cluster ~ "Xa",
                clustered$cluster_id == high_cluster ~ "Xi",
                TRUE ~ NA_character_
            )
        }

        clustered <- clustered %>%
            left_join(cluster_means, by = "cluster_id")
    } else {
        clustered$CLUSTER_MEAN_METHYLATION <- NA_real_
    }

    clustered
}

n_cores <- max(1L, min(threads, parallel::detectCores()))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clustered_list <- foreach(
    i = seq_len(nrow(cpg_bed)),
    .packages = c("NanoMethViz", "dplyr", "tidyr", "readr"),
    .export = c("mbr", "cpg_bed", "min_cluster_reads", "xist_bed", "overlaps_xist")
) %dopar% {
    process_cpg(i)
}

stopCluster(cl)

clustered <- bind_rows(clustered_list)

if (nrow(clustered) == 0) {
    clustered <- tibble(
        read_name = character(),
        cluster_id = character(),
        mean = double(),
        CGI_ID = character(),
        CHROM = character(),
        START = integer(),
        END = integer(),
        assigned_X = character(),
        XIST_LABEL_REVERSED = character(),
        CLUSTER_MEAN_METHYLATION = double()
    )
}

write_tsv(clustered, clustered_out)

assigned <- clustered %>%
    filter(assigned_X %in% c("Xa", "Xi")) %>%
    arrange(CGI_ID, read_name) %>%
    # Follow the published SkewX implementation: a read spanning multiple
    # CpG islands contributes once to the XCI count.
    distinct(read_name, .keep_all = TRUE) %>%
    inner_join(haplotags, by = "read_name") %>%
    filter(HP %in% c(1L, 2L), !is.na(PS), PS != ".")

empty_block_table <- function() {
    tibble(
        PS = character(),
        BLOCK_START = integer(),
        BLOCK_END = integer(),
        CPG_ISLANDS = integer(),
        H1_Xa = integer(),
        H2_Xa = integer(),
        H1_Xi = integer(),
        H2_Xi = integer(),
        TOTAL_READS = integer(),
        H1_XA_SKEW = double(),
        FOLDED_MINOR_XI_FRACTION = double(),
        MAJOR_XI_FRACTION = double(),
        INFORMATIVE_FOR_MLE = character()
    )
}

if (nrow(assigned) == 0) {
    block_skew <- empty_block_table()
} else {
    counts <- assigned %>%
        count(PS, assigned_X, HP, name = "N") %>%
        mutate(
            COMBINATION = paste0(
                "H", HP, "_", assigned_X
            )
        ) %>%
        select(PS, COMBINATION, N) %>%
        pivot_wider(
            names_from = COMBINATION,
            values_from = N,
            values_fill = 0
        )

    for (column in c("H1_Xa", "H2_Xa", "H1_Xi", "H2_Xi")) {
        if (!column %in% names(counts)) {
            counts[[column]] <- 0L
        }
    }

    coordinates <- assigned %>%
        group_by(PS) %>%
        summarise(
            BLOCK_START = min(START, na.rm = TRUE),
            BLOCK_END = max(END, na.rm = TRUE),
            CPG_ISLANDS = n_distinct(CGI_ID),
            .groups = "drop"
        )

    block_skew <- counts %>%
        left_join(coordinates, by = "PS") %>%
        mutate(
            TOTAL_READS = H1_Xa + H2_Xa + H1_Xi + H2_Xi,
            H1_XA_SKEW = ifelse(
                TOTAL_READS > 0,
                (H1_Xa + H2_Xi) / TOTAL_READS,
                NA_real_
            ),
            FOLDED_MINOR_XI_FRACTION = pmin(
                H1_XA_SKEW,
                1 - H1_XA_SKEW
            ),
            MAJOR_XI_FRACTION = 1 - FOLDED_MINOR_XI_FRACTION,
            INFORMATIVE_FOR_MLE = ifelse(
                TOTAL_READS >= min_block_reads,
                "YES",
                "NO"
            )
        ) %>%
        select(
            PS, BLOCK_START, BLOCK_END, CPG_ISLANDS,
            H1_Xa, H2_Xa, H1_Xi, H2_Xi,
            TOTAL_READS, H1_XA_SKEW,
            FOLDED_MINOR_XI_FRACTION,
            MAJOR_XI_FRACTION,
            INFORMATIVE_FOR_MLE
        ) %>%
        arrange(BLOCK_START)
}

write_tsv(block_skew, block_out)

folded_binomial_loglik <- function(p, table) {
    if (p <= 0 || p >= 0.5) {
        # p=0.5 is valid, but keeping the optimisation just inside the
        # boundary avoids numerical edge cases; the reported value may still
        # approach 0.5.
        p <- max(min(p, 0.5 - 1e-10), 1e-10)
    }

    total_ll <- 0
    for (i in seq_len(nrow(table))) {
        n <- table$TOTAL_READS[[i]]
        x <- table$H1_Xa[[i]] + table$H2_Xi[[i]]
        k <- min(x, n - x)

        p1 <- dbinom(k, size = n, prob = p)
        if (k == n - k) {
            probability <- p1
        } else {
            p2 <- dbinom(n - k, size = n, prob = p)
            probability <- p1 + p2
        }
        total_ll <- total_ll + log(max(probability, .Machine$double.xmin))
    }
    total_ll
}

profile_ci <- function(table, mle) {
    if (nrow(table) == 0) {
        return(c(NA_real_, NA_real_))
    }
    grid <- seq(0.001, 0.5, length.out = 2000)
    ll <- vapply(
        grid,
        function(p) folded_binomial_loglik(p, table),
        numeric(1)
    )
    max_ll <- folded_binomial_loglik(mle, table)
    keep <- grid[ll >= (max_ll - 1.920729)]
    if (length(keep) == 0) {
        return(c(NA_real_, NA_real_))
    }
    c(min(keep), max(keep))
}

informative_blocks <- block_skew %>%
    filter(
        INFORMATIVE_FOR_MLE == "YES",
        !is.na(FOLDED_MINOR_XI_FRACTION)
    )

if (nrow(informative_blocks) > 0) {
    fit <- optimize(
        function(p) -folded_binomial_loglik(p, informative_blocks),
        interval = c(1e-6, 0.5),
        tol = 1e-7
    )
    mle_minor <- fit$minimum
    ci <- profile_ci(informative_blocks, mle_minor)
    mle_major <- 1 - mle_minor
    ratio <- sprintf("%.3f:%.3f", mle_minor, mle_major)
    status <- "XCI_SKEW_ESTIMATED"
} else {
    mle_minor <- NA_real_
    mle_major <- NA_real_
    ci <- c(NA_real_, NA_real_)
    ratio <- "."
    status <- "INSUFFICIENT_INFORMATIVE_PHASE_BLOCKS"
}

summary <- tibble(
    SAMPLE = sample_id,
    XCI_STATUS = status,
    CPG_ISLANDS_INPUT = nrow(cpg_bed),
    CPG_ISLANDS_WITH_TWO_METHYLATION_CLUSTERS = n_distinct(
        clustered$CGI_ID[clustered$assigned_X %in% c("Xa", "Xi")]
    ),
    PHASED_METHYLATION_READS = nrow(assigned),
    INFORMATIVE_PHASE_BLOCKS = nrow(informative_blocks),
    INFORMATIVE_READS = sum(informative_blocks$TOTAL_READS, na.rm = TRUE),
    XCI_MINOR_INACTIVE_FRACTION_MLE = mle_minor,
    XCI_MAJOR_INACTIVE_FRACTION_MLE = mle_major,
    XCI_MINOR_FRACTION_CI95_LOW = ci[[1]],
    XCI_MINOR_FRACTION_CI95_HIGH = ci[[2]],
    XCI_FOLDED_RATIO = ratio,
    HIGH_SKEW_80_20_OR_MORE = ifelse(
        !is.na(mle_minor) && mle_minor <= 0.20,
        "YES",
        "NO"
    ),
    EXTREME_SKEW_90_10_OR_MORE = ifelse(
        !is.na(mle_minor) && mle_minor <= 0.10,
        "YES",
        "NO"
    ),
    XIST_PROMOTER_SPECIAL_HANDLING = ifelse(
        is.null(xist_bed),
        "NOT_CONFIGURED",
        "CONFIGURED_REVERSED_XA_XI_LABEL"
    ),
    XCI_METHOD = (
        "CpG-island read methylation clustering + WhatsHap HP/PS + "
        "folded-binomial maximum-likelihood skew estimate"
    ),
    INTERPRETATION = (
        "0.5 indicates balanced XCI; smaller minor inactive-X fractions "
        "indicate stronger skew. Haplotype identity can flip between phase "
        "blocks, so this folded estimate measures degree of skew, not parental "
        "directionality."
    )
)

write_tsv(summary, summary_out)

message(
    "[OK] XCI sample=", sample_id,
    " status=", status,
    " blocks=", nrow(informative_blocks),
    " minor_fraction=", ifelse(is.na(mle_minor), ".", sprintf("%.4f", mle_minor))
)
