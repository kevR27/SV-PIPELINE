library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
library(svglite)

snakemake@source("theme_thesis.R")

bins_file <- snakemake@input[["bins"]]
summary_file <- snakemake@input[["summary"]]
genes_file <- snakemake@input[["genes"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_n <- as.integer(snakemake@params[["top_n"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

bins <- fread(bins_file, na.strings = c("", ".", "NA"))
summary <- fread(summary_file, na.strings = c("", ".", "NA"))
genes <- fread(genes_file, na.strings = c("", ".", "NA"))

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

summary[, SIZE_MB := as.numeric(SV_SIZE_BP) / 1e6]
summary <- head(summary[order(-SIZE_MB)], top_n)

lead_col <- if ("LEAD_GENE" %in% names(summary)) "LEAD_GENE" else "GENE"
summary[, LABEL := paste0(
  get(lead_col), " | ", SVTYPE, " | ", sprintf("%.2f Mb", SIZE_MB), "\n",
  "depth ratio=", sprintf("%.2f", as.numeric(DEPTH_RATIO)),
  " | ", DEPTH_PATTERN
)]

bins <- bins[SV_ID %in% summary$SV_ID]
bins <- merge(
  bins,
  summary[, .(SV_ID, DEPTH_RATIO, DEPTH_PATTERN, SIZE_MB, LABEL)],
  by = "SV_ID",
  all.x = TRUE
)
bins[, POSITION_MB := (as.numeric(PLOT_START) + as.numeric(PLOT_END)) / 2 / 1e6]
setorder(bins, SV_ID, PLOT_START)
bins[, ROLLING_MEDIAN_DEPTH := frollmedian(
  as.numeric(NORMALIZED_DEPTH),
  n = 3,
  align = "center",
  fill = NA_real_,
  na.rm = TRUE
), by = SV_ID]

boundary <- unique(bins[, .(SV_ID, LABEL, SV_START, SV_END)])
boundary[, SV_START_MB := as.numeric(SV_START) / 1e6]
boundary[, SV_END_MB := as.numeric(SV_END) / 1e6]

gene_plot <- genes[SV_ID %in% summary$SV_ID & !is.na(GENE_MID)]
if (nrow(gene_plot) > 0) {
  gene_plot <- merge(
    gene_plot,
    summary[, .(SV_ID, LABEL)],
    by = "SV_ID",
    all.x = TRUE
  )
  gene_plot[, GENE_MID_MB := as.numeric(GENE_MID) / 1e6]
  gene_plot[, SCORE_NUM := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
  gene_plot[is.na(SCORE_NUM), SCORE_NUM := 0]
  gene_plot[, PANEL_FLAG := fifelse(PANEL_STATUS == "PANEL_GENE", 1L, 0L)]
  gene_plot[, MITO_FLAG := fifelse(
    MITOCARTA %in% c("NUCLEAR_MITOCHONDRIAL_GENE", "MTDNA_ENCODED_GENE"),
    1L,
    0L
  )]
  gene_plot[, LABEL_PRIORITY := PANEL_FLAG * 100 + MITO_FLAG * 10 + SCORE_NUM]

  label_genes <- gene_plot[
    order(-LABEL_PRIORITY, GENE),
    head(.SD, 8),
    by = SV_ID
  ]
} else {
  label_genes <- data.table()
}

p <- ggplot(bins, aes(POSITION_MB, as.numeric(NORMALIZED_DEPTH))) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_line(linewidth = 0.35, colour = "grey70", alpha = 0.6) +
  geom_line(
    aes(y = ROLLING_MEDIAN_DEPTH),
    linewidth = 0.9,
    colour = "grey15",
    na.rm = TRUE
  ) +
  geom_vline(
    data = boundary,
    aes(xintercept = SV_START_MB),
    inherit.aes = FALSE,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey45"
  ) +
  geom_vline(
    data = boundary,
    aes(xintercept = SV_END_MB),
    inherit.aes = FALSE,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey45"
  )

if (nrow(gene_plot) > 0) {
  p <- p +
    geom_segment(
      data = gene_plot,
      aes(
        x = GENE_MID_MB,
        xend = GENE_MID_MB,
        y = -0.16,
        yend = -0.04
      ),
      inherit.aes = FALSE,
      linewidth = 0.28,
      colour = "grey50",
      alpha = 0.75
    )
}

if (nrow(label_genes) > 0) {
  p <- p +
    geom_text_repel(
      data = label_genes,
      aes(x = GENE_MID_MB, y = -0.18, label = GENE),
      inherit.aes = FALSE,
      size = 2.5,
      direction = "x",
      angle = 35,
      hjust = 0,
      segment.size = 0.2,
      min.segment.length = 0,
      box.padding = 0.18,
      point.padding = 0.05,
      max.overlaps = Inf,
      seed = 1
    )
}

p <- p +
  facet_wrap(~ LABEL, scales = "free_x", ncol = 2) +
  scale_y_continuous(
    name = "Normalized binned depth",
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.13, 0.08))
  ) +
  scale_x_continuous(
    name = "Genomic position (Mb)",
    labels = label_number(accuracy = 0.1)
  ) +
  expand_limits(y = -0.28) +
  labs(
    title = "Read-depth support and gene context for large copy-number-changing SVs",
    subtitle = paste(
      "Grey: adaptive median depth bins; black: 3-bin rolling median; dashed lines: SV boundaries.",
      "Gene ticks are aligned to genomic position; labels prioritize ON-panel, MitoCarta and high-relevance genes."
    )
  ) +
  theme_thesis(10) +
  theme(
    strip.text = element_text(size = 8.2, lineheight = 0.95),
    panel.spacing = unit(1.0, "lines"),
    plot.subtitle = element_text(size = 9.2, margin = margin(b = 8))
  )

height <- max(7, ceiling(nrow(summary) / 2) * 3.1 + 1.8)
ggsave(pdf_file, p, width = 14.5, height = height, device = cairo_pdf)
ggsave(png_file, p, width = 14.5, height = height, dpi = 400)
ggsave(svg_file, p, width = 14.5, height = height, device = svglite)
