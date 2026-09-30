library(data.table)
library(ggplot2)
library(scales)
library(svglite)

snakemake@source("theme_thesis.R")

gene_file <- snakemake@input[["genes"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_n <- as.integer(snakemake@params[["top_n"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

dt <- fread(gene_file, na.strings = c("", ".", "NA", "N/A"))

if (nrow(dt) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No population-evidence gene rows available", size = 5) +
    theme_void()
} else {
  dt[, SCORE := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
  dt[is.na(SCORE), SCORE := 0]
  setorder(dt, -SCORE, -RARE_OR_NOT_OBSERVED_NEEDLR_SV_COUNT, COMMON_NEEDLR_SV_COUNT, GENE)
  dt <- head(dt, top_n)

  long <- melt(
    dt,
    id.vars = c("GENE", "PANEL_STATUS", "SCORE"),
    measure.vars = c(
      "RARE_OR_NOT_OBSERVED_NEEDLR_SV_COUNT",
      "COMMON_NEEDLR_SV_COUNT",
      "GNOMAD_BENIGN_OVERLAP_SV_COUNT"
    ),
    variable.name = "EVIDENCE",
    value.name = "SV_COUNT"
  )

  labels <- c(
    RARE_OR_NOT_OBSERVED_NEEDLR_SV_COUNT = "Rare/not observed in needLR",
    COMMON_NEEDLR_SV_COUNT = "Common in needLR",
    GNOMAD_BENIGN_OVERLAP_SV_COUNT = "gnomAD included in AnnotSV benign overlap"
  )
  long[, EVIDENCE := labels[as.character(EVIDENCE)]]
  long[, GENE := factor(GENE, levels = rev(unique(dt$GENE)))]

  p <- ggplot(long, aes(x = GENE, y = SV_COUNT, fill = EVIDENCE)) +
    geom_col(position = "dodge", width = 0.78) +
    coord_flip() +
    labs(
      title = "Population-frequency context of prioritized SV-overlapping genes",
      subtitle = paste(
        "needLR represents matched ONT control-frequency evidence.",
        "The gnomAD bar indicates presence within AnnotSV benign-region sources; it is not a gnomAD-specific exact-allele AF."
      ),
      x = NULL,
      y = "Unique structural variants",
      fill = NULL
    ) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_thesis(10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.subtitle = element_text(size = 9, lineheight = 1.0, margin = margin(b = 8)),
      axis.text.y = element_text(size = 8.5)
    )
}

height <- max(7.5, min(13, top_n * 0.34 + 2.5))
ggsave(pdf_file, p, width = 12.5, height = height, device = cairo_pdf)
ggsave(png_file, p, width = 12.5, height = height, dpi = 400)
ggsave(svg_file, p, width = 12.5, height = height, device = svglite)
