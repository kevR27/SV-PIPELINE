library(data.table)
library(ggplot2)
library(scales)
library(svglite)

snakemake@source("theme_thesis.R")

bins_file <- snakemake@input[["bins"]]
summary_file <- snakemake@input[["summary"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_n <- as.integer(snakemake@params[["top_n"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

bins <- fread(bins_file, na.strings = c("", ".", "NA"))
summary <- fread(summary_file, na.strings = c("", ".", "NA"))
if (nrow(bins) == 0 || nrow(summary) == 0) {
  empty_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No DEL/DUP >= configured size threshold", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = "Read-depth support for large copy-number-changing SVs") +
    theme_void()
  ggsave(pdf_file, empty_plot, width = 10, height = 4.5)
  ggsave(png_file, empty_plot, width = 10, height = 4.5, dpi = 400)
  ggsave(svg_file, empty_plot, width = 10, height = 4.5, device = svglite)
  quit(save = "no", status = 0)
}

summary[, SIZE_MB := SV_SIZE_BP / 1e6]
summary <- head(summary[order(-SIZE_MB)], top_n)
bins <- bins[SV_ID %in% summary$SV_ID]
bins <- merge(bins, summary[, .(SV_ID, DEPTH_RATIO, DEPTH_PATTERN, SIZE_MB)], by = "SV_ID", all.x = TRUE)
bins[, POSITION_MB := (PLOT_START + PLOT_END) / 2 / 1e6]
setorder(bins, SV_ID, PLOT_START)
bins[, ROLLING_MEDIAN_DEPTH := frollmedian(
  NORMALIZED_DEPTH,
  n = 3,
  align = "center",
  fill = NA_real_,
  na.rm = TRUE
), by = SV_ID]

bins[, LABEL := paste0(
  GENE, " | ", SVTYPE, " | ", sprintf("%.2f Mb", SIZE_MB),
  " | depth ratio=", sprintf("%.2f", as.numeric(DEPTH_RATIO))
)]

boundary <- unique(bins[, .(LABEL, SV_START, SV_END)])
boundary[, SV_START_MB := SV_START / 1e6]
boundary[, SV_END_MB := SV_END / 1e6]

p <- ggplot(bins, aes(POSITION_MB, NORMALIZED_DEPTH)) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.35, colour = "grey65", alpha = 0.55) +
  geom_line(
    aes(y = ROLLING_MEDIAN_DEPTH),
    linewidth = 0.85,
    colour = "grey15",
    na.rm = TRUE
  ) +
  geom_vline(data = boundary, aes(xintercept = SV_START_MB), linewidth = 0.4, linetype = "dashed", colour = "grey45") +
  geom_vline(data = boundary, aes(xintercept = SV_END_MB), linewidth = 0.4, linetype = "dashed", colour = "grey45") +
  facet_wrap(~ LABEL, scales = "free_x", ncol = 2) +
  scale_y_continuous(
    name = "Normalized binned depth",
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  scale_x_continuous(
    name = "Genomic position (Mb)",
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Read-depth support for large copy-number-changing SVs",
    subtitle = "Grey: adaptive median depth bins. Black: 3-bin rolling median. Dashed lines mark SV boundaries."
  ) +
  theme_thesis(10) +
  theme(strip.text = element_text(size = 8.5))

height <- max(7, ceiling(nrow(summary) / 2) * 2.5 + 1.5)
ggsave(pdf_file, p, width = 13, height = height, device = cairo_pdf)
ggsave(png_file, p, width = 13, height = height, dpi = 400)
ggsave(svg_file, p, width = 13, height = height, device = svglite)
