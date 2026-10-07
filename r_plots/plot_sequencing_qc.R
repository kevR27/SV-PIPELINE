#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

qc_file <- snakemake@input[["qc"]]
n50_file <- snakemake@input[["n50"]]

pdf_out <- snakemake@output[["pdf"]]
png_out <- snakemake@output[["png"]]
svg_out <- snakemake@output[["svg"]]

qc <- fread(qc_file, na.strings = c("", ".", "NA"))
n50 <- fread(n50_file, na.strings = c("", ".", "NA"))

if (nrow(qc) == 0) {
  stop("Dorado QC table is empty")
}

sample_id <- as.character(qc$SAMPLE[[1]])

theme_thesis <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  )

# Read-length statistics from the N50 table.
length_metrics <- melt(
  n50[
    metric %in% c("mean_read_length", "median_read_length", "N50")
  ],
  id.vars = "metric",
  variable.name = "read_group",
  value.name = "value"
)
length_metrics[, value := as.numeric(value)]
length_metrics[, metric := factor(
  metric,
  levels = c("mean_read_length", "median_read_length", "N50"),
  labels = c("Mean", "Median", "N50")
)]

p_length <- ggplot(
  length_metrics[is.finite(value)],
  aes(read_group, value, fill = metric)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = " kb")) +
  labs(
    title = "Read length",
    x = NULL,
    y = "Length",
    fill = NULL
  ) +
  theme_thesis +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "top"
  )

# Long-read vs short-read contribution.
fraction_metrics <- melt(
  n50[
    metric %in% c("read_fraction_percent", "base_fraction_percent")
  ],
  id.vars = "metric",
  variable.name = "read_group",
  value.name = "percent"
)
fraction_metrics[, percent := as.numeric(percent)]
fraction_metrics <- fraction_metrics[
  read_group != "all_primary_reads"
]
fraction_metrics[, metric := factor(
  metric,
  levels = c("read_fraction_percent", "base_fraction_percent"),
  labels = c("Reads", "Bases")
)]

p_fraction <- ggplot(
  fraction_metrics[is.finite(percent)],
  aes(read_group, percent, fill = metric)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Read and base contribution",
    x = NULL,
    y = "Percent",
    fill = NULL
  ) +
  theme_thesis +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "top"
  )

# Q-score summary.
qscore <- data.table(
  metric = c("Mean Q-score", "Median Q-score"),
  value = c(
    as.numeric(qc$MEAN_QSCORE[[1]]),
    as.numeric(qc$MEDIAN_QSCORE[[1]])
  )
)

p_qscore <- ggplot(qscore[is.finite(value)], aes(metric, value)) +
  geom_col(width = 0.55) +
  labs(
    title = "Basecalling quality",
    x = NULL,
    y = "Q-score"
  ) +
  theme_thesis

# Main sequencing metrics as a compact text panel.
fmt_int <- function(x) {
  ifelse(is.na(x), ".", format(round(as.numeric(x)), big.mark = ",", scientific = FALSE))
}
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), ".", format(round(as.numeric(x), digits), nsmall = digits))
}

passing <- if ("PASSING_READ_PERCENT" %in% names(qc)) {
  paste0(fmt_num(qc$PASSING_READ_PERCENT[[1]], 1), "%")
} else {
  "."
}

summary_lines <- c(
  paste0("Reads: ", fmt_int(qc$READ_COUNT[[1]])),
  paste0("Yield: ", fmt_num(qc$YIELD_GBP[[1]], 2), " Gb"),
  paste0("Read N50: ", fmt_num(qc$READ_N50[[1]] / 1000, 2), " kb"),
  paste0("Mean read length: ", fmt_num(qc$MEAN_READ_LENGTH[[1]] / 1000, 2), " kb"),
  paste0("Median read length: ", fmt_num(qc$MEDIAN_READ_LENGTH[[1]] / 1000, 2), " kb"),
  paste0("Mean Q-score: ", fmt_num(qc$MEAN_QSCORE[[1]], 2)),
  paste0("Passing reads: ", passing)
)

p_summary <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = seq(length(summary_lines), 1),
    label = summary_lines,
    hjust = 0,
    size = 4
  ) +
  xlim(0, 1) +
  ylim(0.5, length(summary_lines) + 0.5) +
  labs(title = "Sequencing summary") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  )

combined <- (
  p_summary | p_qscore
) / (
  p_length | p_fraction
) +
  plot_annotation(
    title = paste0(sample_id, ": Oxford Nanopore sequencing QC"),
    subtitle = "Dorado summary and BAM read-length metrics"
  )

dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)
ggsave(pdf_out, combined, width = 13.5, height = 9.0, units = "in")
ggsave(png_out, combined, width = 13.5, height = 9.0, units = "in", dpi = 300)
ggsave(svg_out, combined, width = 13.5, height = 9.0, units = "in")

message("[OK] sequencing QC plot written for ", sample_id)
