#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(
    paste(
      "Usage: plot_xci.R",
      "<block_skew.tsv> <summary.tsv> <phase_blocks.tsv> <phase_summary.tsv>",
      "<hp1.bedmethyl.gz> <hp2.bedmethyl.gz> <out_dir> <sample> <bin_bp>"
    )
  )
}

block_file <- args[[1]]
summary_file <- args[[2]]
phase_file <- args[[3]]
phase_summary_file <- args[[4]]
hp1_file <- args[[5]]
hp2_file <- args[[6]]
out_dir <- args[[7]]
sample_id <- args[[8]]
bin_bp <- as.numeric(args[[9]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, width, height) {
  prefix <- file.path(out_dir, paste0(sample_id, "_", name))
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, units = "in")
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height, units = "in")
}

theme_thesis <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10.5),
      legend.position = "top"
    )
}

blocks <- fread(block_file, na.strings = c("", ".", "NA"))
summary <- fread(summary_file, na.strings = c("", ".", "NA"))
phase <- fread(phase_file, na.strings = c("", ".", "NA"))
phase_summary <- fread(phase_summary_file, na.strings = c("", ".", "NA"))

global_p <- if (nrow(summary) > 0 && "GLOBAL_FOLDED_SKEW_P" %in% names(summary)) {
  as.numeric(summary$GLOBAL_FOLDED_SKEW_P[[1]])
} else {
  NA_real_
}

ratio <- if (nrow(summary) > 0 && "XCI_MAJOR_MINOR_RATIO" %in% names(summary)) {
  as.character(summary$XCI_MAJOR_MINOR_RATIO[[1]])
} else {
  "."
}

status <- if (nrow(summary) > 0 && "XCI_SKEW_STATUS" %in% names(summary)) {
  as.character(summary$XCI_SKEW_STATUS[[1]])
} else {
  "."
}

# 1. Global skew distribution.
if (nrow(blocks) > 0) {
  blocks[, FOLDED_BLOCK_SKEW := as.numeric(FOLDED_BLOCK_SKEW)]
  blocks[, TRIALS := as.numeric(TRIALS)]

  p1 <- ggplot(blocks[is.finite(FOLDED_BLOCK_SKEW) & TRIALS > 0]) +
    geom_histogram(
      aes(x = FOLDED_BLOCK_SKEW, weight = TRIALS),
      binwidth = 0.025,
      boundary = 0
    ) +
    geom_vline(xintercept = 0.5, linetype = "dotted") +
    labs(
      title = paste0(sample_id, ": X-inactivation skew"),
      subtitle = paste0("Global estimate: ", ratio, " | ", status),
      x = "Folded X-inactivation skew",
      y = "Supporting reads"
    ) +
    coord_cartesian(xlim = c(0, 0.5)) +
    theme_thesis()

  if (is.finite(global_p)) {
    p1 <- p1 + geom_vline(xintercept = global_p, linetype = "dashed", linewidth = 0.9)
  }
} else {
  p1 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No informative X-inactivation blocks") +
    theme_void() +
    labs(title = paste0(sample_id, ": X-inactivation skew"))
}
save_plot(p1, "xci_block_skew_distribution", 9.5, 6.0)

# 2. Local skew along chrX.
if (nrow(blocks) > 0) {
  blocks[, BLOCK_START := as.numeric(BLOCK_START)]
  blocks[, BLOCK_END := as.numeric(BLOCK_END)]
  blocks[, H1_Xa_SKEW := as.numeric(H1_Xa_SKEW)]
  blocks[, BLOCK_MID_MB := (BLOCK_START + BLOCK_END) / 2 / 1e6]

  p2 <- ggplot(blocks[is.finite(H1_Xa_SKEW)], aes(BLOCK_MID_MB, H1_Xa_SKEW)) +
    geom_point(aes(size = TRIALS), alpha = 0.75) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    scale_size_continuous(name = "Informative reads", range = c(2, 8)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(sample_id, ": X-inactivation along chrX"),
      subtitle = "Values above 0.5 support H1 as Xa; values below 0.5 support H2 as Xa",
      x = "chrX position (Mb)",
      y = "H1 Xa proportion"
    ) +
    theme_thesis()
} else {
  p2 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No informative X-inactivation blocks") +
    theme_void() +
    labs(title = paste0(sample_id, ": X-inactivation along chrX"))
}
save_plot(p2, "xci_chrX_block_skew", 12.0, 5.8)

# 3. Confidence in Xa/Xi orientation.
if (
  nrow(blocks) > 0 &&
  "LOG10_ODDS_H1_XA_VS_H2_XA" %in% names(blocks)
) {
  blocks[, LOG10_ODDS_H1_XA_VS_H2_XA := as.numeric(LOG10_ODDS_H1_XA_VS_H2_XA)]

  p3 <- ggplot(
    blocks[is.finite(LOG10_ODDS_H1_XA_VS_H2_XA) & TRIALS > 0],
    aes(TRIALS, LOG10_ODDS_H1_XA_VS_H2_XA)
  ) +
    geom_point(aes(size = TRIALS), alpha = 0.75) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted") +
    scale_size_continuous(name = "Informative reads", range = c(2, 8)) +
    labs(
      title = paste0(sample_id, ": confidence in X-inactivation orientation"),
      subtitle = "|log10 odds| >= 1 means at least 10:1 support; >= 2 means at least 100:1 support",
      x = "Informative reads in phase block",
      y = "log10 odds: H1 Xa vs H2 Xa"
    ) +
    theme_thesis()
} else {
  p3 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No block-level orientation odds available") +
    theme_void() +
    labs(title = paste0(sample_id, ": confidence in X-inactivation orientation"))
}
save_plot(p3, "xci_orientation_log_odds", 10.5, 6.0)

read_bedmethyl <- function(path, haplotype) {
  cmd <- paste("gzip -dc", shQuote(path))
  dt <- fread(cmd = cmd, header = FALSE, select = 1:11, showProgress = FALSE)
  if (ncol(dt) < 11) {
    stop(paste("bedMethyl has fewer than 11 columns:", path))
  }
  setnames(
    dt,
    c(
      "chrom", "start", "end", "name", "score", "strand",
      "thick_start", "thick_end", "color", "coverage", "methylation"
    )
  )
  dt[, start := as.numeric(start)]
  dt[, coverage := as.numeric(coverage)]
  dt[, methylation := as.numeric(methylation)]
  dt <- dt[
    is.finite(start) &
    is.finite(coverage) &
    is.finite(methylation)
  ]
  dt[, bin_start := floor(start / bin_bp) * bin_bp]
  dt[, weighted_methylation := methylation * coverage]

  out <- dt[
    ,
    .(
      coverage_sum = sum(coverage, na.rm = TRUE),
      weighted_sum = sum(weighted_methylation, na.rm = TRUE)
    ),
    by = bin_start
  ]
  out[, methylation_percent := weighted_sum / coverage_sum]
  out[, position_mb := (bin_start + bin_bp / 2) / 1e6]
  out[, haplotype := haplotype]
  out
}

hp1 <- read_bedmethyl(hp1_file, "HP1")
hp2 <- read_bedmethyl(hp2_file, "HP2")
meth <- rbindlist(list(hp1, hp2), fill = TRUE)

# 4. Haplotype-specific chrX methylation.
if (nrow(meth) > 0) {
  p4 <- ggplot(
    meth[is.finite(methylation_percent)],
    aes(position_mb, methylation_percent, group = haplotype, linetype = haplotype)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      title = paste0(sample_id, ": phased X-chromosome methylation"),
      subtitle = "CpG 5mC summarized separately for WhatsHap HP1 and HP2",
      x = "chrX position (Mb)",
      y = "CpG 5mC (%)",
      linetype = "Haplotype"
    ) +
    theme_thesis()
} else {
  p4 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No haplotype-specific chrX methylation available") +
    theme_void() +
    labs(title = paste0(sample_id, ": phased X-chromosome methylation"))
}
save_plot(p4, "xci_haplotype_methylation", 12.0, 5.8)

# 5. WhatsHap vs LongPhase phase consistency.
if (nrow(phase) > 0) {
  phase[, ORIENTATION_CONCORDANCE := as.numeric(ORIENTATION_CONCORDANCE)]
  phase[, SHARED_PHASED_SNVS := as.numeric(SHARED_PHASED_SNVS)]
  phase <- phase[order(-SHARED_PHASED_SNVS, -ORIENTATION_CONCORDANCE)]
  phase <- head(phase, 30)
  phase[, block := paste0("W", WHATSHAP_PS, "/L", LONGPHASE_PS)]
  phase[, block := factor(block, levels = rev(block))]

  overall <- if (
    nrow(phase_summary) > 0 &&
    "FLIP_TOLERANT_PHASE_CONCORDANCE" %in% names(phase_summary)
  ) {
    as.character(phase_summary$FLIP_TOLERANT_PHASE_CONCORDANCE[[1]])
  } else {
    "."
  }

  p5 <- ggplot(phase, aes(block, ORIENTATION_CONCORDANCE)) +
    geom_col() +
    coord_flip(ylim = c(0.5, 1)) +
    labs(
      title = paste0(sample_id, ": WhatsHap and LongPhase agreement"),
      subtitle = paste0("Overall flip-tolerant concordance: ", overall),
      x = "Overlapping phase blocks",
      y = "Phase concordance"
    ) +
    theme_thesis()
} else {
  p5 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No shared chrX phase blocks") +
    theme_void() +
    labs(title = paste0(sample_id, ": WhatsHap and LongPhase agreement"))
}
save_plot(p5, "xci_phase_concordance", 10.5, 6.5)

manifest <- data.table(
  sample = sample_id,
  XCI_SKEW_STATUS = status,
  XCI_MAJOR_MINOR_RATIO = ratio,
  GLOBAL_FOLDED_SKEW_P = ifelse(is.finite(global_p), global_p, NA_real_),
  plot_directory = out_dir
)
fwrite(
  manifest,
  file.path(out_dir, paste0(sample_id, "_xci_plot_manifest.tsv")),
  sep = "\t",
  na = "."
)

message("[OK] X-inactivation R plots written to ", out_dir)
