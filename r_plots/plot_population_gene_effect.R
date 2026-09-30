library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(svglite)

snakemake@source("theme_thesis.R")

input_file <- snakemake@input[["table"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_n <- as.integer(snakemake@params[["top_n"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file, na.strings = c("", ".", "NA"))
if (nrow(dt) == 0) stop("No population/effect rows available.")

dt <- unique(dt, by = c("SV_ID", "GENE"))
dt[, SCORE_NUM := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
dt[is.na(SCORE_NUM), SCORE_NUM := 0]
dt[, RARE_FLAG := as.integer(POPULATION_CLASS == "LOW_FREQUENCY_BY_NEEDLR")]
dt[, PANEL_FLAG := as.integer(PANEL_STATUS == "PANEL_GENE")]

gene_rank <- dt[, .(
  RARE_COUNT = sum(RARE_FLAG),
  PANEL_FLAG = max(PANEL_FLAG),
  SCORE = max(SCORE_NUM),
  SV_COUNT = uniqueN(SV_ID)
), by = GENE]
setorder(gene_rank, -PANEL_FLAG, -RARE_COUNT, -SCORE, -SV_COUNT, GENE)
selected <- head(gene_rank$GENE, top_n)
plot_dt <- dt[GENE %in% selected]
plot_dt[, GENE := factor(GENE, levels = rev(selected))]

population_labels <- c(
  "LOW_FREQUENCY_BY_NEEDLR" = "Low frequency\nneedLR",
  "COMMON_BY_NEEDLR" = "Common\nneedLR",
  "COMMON_BENIGN_REGION_OVERLAP_CONTEXT" = "Common benign-region\noverlap",
  "LOW_AF_BENIGN_REGION_OVERLAP_CONTEXT" = "Low-AF benign-region\noverlap",
  "NO_NEEDLR_MATCH_AF_UNKNOWN" = "No needLR match\nAF unknown",
  "POPULATION_AF_NOT_EVALUABLE" = "AF not\nevaluable",
  "POPULATION_FREQUENCY_UNKNOWN" = "Frequency\nunknown"
)
effect_labels <- c(
  "DIRECT_BREAKPOINT" = "Direct\nbreakpoint",
  "COPY_LOSS_GEOMETRY" = "Copy-loss\ngeometry",
  "COPY_GAIN_GEOMETRY" = "Copy-gain\ngeometry",
  "INSERTION_IN_GENE" = "Insertion\nin gene",
  "INVERSION_SPANNED_GENE" = "Inversion-spanned\ngene",
  "NEAR_GENE" = "Near\ngene",
  "INVERSION_OTHER" = "Other\ninversion",
  "OTHER_OR_UNRESOLVED" = "Other /\nunresolved"
)

plot_dt[, POP_LABEL := fifelse(
  POPULATION_CLASS %in% names(population_labels),
  population_labels[POPULATION_CLASS],
  gsub("_", " ", POPULATION_CLASS)
)]
plot_dt[, EFFECT_LABEL := fifelse(
  SV_EFFECT_GROUP %in% names(effect_labels),
  effect_labels[SV_EFFECT_GROUP],
  gsub("_", " ", SV_EFFECT_GROUP)
)]

pop <- plot_dt[, .(N = uniqueN(SV_ID)), by = .(GENE, POP_LABEL)]
effect <- plot_dt[, .(N = uniqueN(SV_ID)), by = .(GENE, EFFECT_LABEL)]

p1 <- ggplot(pop, aes(POP_LABEL, GENE, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(N > 0, N, "")), size = 2.5) +
  scale_fill_gradient(low = "grey95", high = "grey20") +
  labs(
    title = "Population-frequency context",
    x = NULL, y = NULL, fill = "SV count"
  ) +
  theme_thesis(9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right"
  )

p2 <- ggplot(effect, aes(EFFECT_LABEL, GENE, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(N > 0, N, "")), size = 2.5) +
  scale_fill_gradient(low = "grey95", high = "grey20") +
  labs(
    title = "Predicted SV-gene relationship",
    x = NULL, y = NULL, fill = "SV count"
  ) +
  theme_thesis(9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

combined <- p1 + p2 +
  plot_annotation(
    title = "Population context and structural-variant effects by gene",
    subtitle = paste(
      "needLR provides matched long-read control-frequency evidence. gnomAD-SV is shown as AnnotSV benign-source overlap;",
      "AnnotSV benign AFmax may combine databases and is not treated as a gnomAD-specific exact-allele AF."
    )
  )

height <- max(8, 0.28 * length(selected) + 3.2)
ggsave(pdf_file, combined, width = 16, height = height, device = cairo_pdf)
ggsave(png_file, combined, width = 16, height = height, dpi = 400)
ggsave(svg_file, combined, width = 16, height = height, device = svglite)
