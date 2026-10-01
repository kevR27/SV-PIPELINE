library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)

input_file <- snakemake@input[["candidates"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_n <- as.integer(snakemake@params[["top_n"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file, na.strings = c("", ".", "NA", "N/A"))
if (nrow(dt) == 0) stop("No SV-gene candidates available for heatmap.")

dt[, CALLER_COUNT_NUM := suppressWarnings(as.numeric(CALLER_COUNT))]
dt[is.na(CALLER_COUNT_NUM), CALLER_COUNT_NUM := 0]
dt[, SCORE_NUM := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
dt[is.na(SCORE_NUM), SCORE_NUM := 0]
setorder(dt, -SCORE_NUM, -CALLER_COUNT_NUM)
dt <- unique(dt, by = c("SV_ID", "GENE"))
dt <- head(dt, top_n)

flag <- function(x, positive, unknown = character()) {
  x <- as.character(x)
  ifelse(
    is.na(x) | x %in% unknown,
    NA_real_,
    ifelse(x %in% positive, 1, 0)
  )
}

functional_context <- if ("SV_FUNCTIONAL_CONTEXT" %in% names(dt)) {
  dt$SV_FUNCTIONAL_CONTEXT
} else {
  rep("", nrow(dt))
}
effect_direct <- grepl(
  "DIRECT_TRANSCRIPT_DISRUPTION|COPY_LOSS_GEOMETRIC_CONTEXT|COPY_GAIN_GEOMETRIC_CONTEXT|INSERTION_SITE_CONTEXT",
  functional_context
)
effect_regulatory <- grepl(
  "REGULATORY_OR_POSITION_EFFECT|FULLY_SPANNED_COPY_NEUTRAL_REGULATORY_3D",
  functional_context
)
depth_support <- if ("DEPTH_SUPPORT_CLASS" %in% names(dt)) {
  flag(
    dt$DEPTH_SUPPORT_CLASS,
    "SUPPORTS_CALLED_COPY_CHANGE",
    c("NOT_APPLICABLE", "NO_DEPTH_RESULT", "NO_DEPTH_SUMMARY", "REVIEW")
  )
} else {
  rep(NA_real_, nrow(dt))
}

evidence <- cbind(
  flag(dt$CALL_SUPPORT, "MULTI_CALLER", c("UNKNOWN")),
  flag(
    dt$POPULATION_STATUS,
    "PROVISIONAL_LOW_FREQUENCY",
    c("UNKNOWN", "NO_POPULATION_MATCH", "NOT_EVALUABLE_GE_10MB", "NOT_EVALUABLE_BREAKEND")
  ),
  flag(dt$PANEL_STATUS, "PANEL_GENE", c("UNKNOWN")),
  flag(dt$MITOCARTA, c("NUCLEAR_MITOCHONDRIAL_GENE", "MTDNA_ENCODED_GENE"), c("UNKNOWN")),
  ifelse(effect_direct, 1, 0),
  ifelse(effect_regulatory, 1, 0),
  depth_support,
  flag(dt$LONGPHASE_PHASED, "YES", c("NOT_AVAILABLE", "UNKNOWN")),
  flag(dt$STRAGLR, "YES", c("NOT_AVAILABLE", "UNKNOWN")),
  flag(dt$TLDR, "YES", c("NOT_AVAILABLE", "UNKNOWN"))
)
colnames(evidence) <- c(
  "Multi-caller", "Low-frequency needLR match", "ON panel", "MitoCarta",
  "Direct/geometry SV-gene relation", "Regulatory/inversion context",
  "Depth supports copy change", "LongPhase phased",
  "Straglr same-locus insertion", "TLDR insertion match"
)

rownames(evidence) <- paste0(
  dt$GENE, " | ", dt$SVTYPE, " | ", dt$CHROM, ":",
  format(as.numeric(dt$START), scientific = FALSE, trim = TRUE)
)
group_labels <- c(
  "COPY_NUMBER_SV_GE_10MB" = "Copy-number-changing SV >=10 Mb",
  "LARGE_COPY_NUMBER_SV" = "Large copy-number-changing SV",
  "COPY_NUMBER_SV" = "Copy-number-changing SV",
  "INSERTION" = "Insertion",
  "BREAKPOINT_IN_GENE" = "Breakpoint in gene",
  "BREAKPOINT_NEAR_GENE" = "Breakpoint near gene",
  "INVERSION_SPANNED_GENE" = "Inversion-spanned gene",
  "GENE_INSIDE_REARRANGEMENT" = "Inversion-spanned gene",
  "BREAKPOINT_SV" = "Breakpoint SV",
  "OTHER_SV" = "Other SV"
)
group_text <- as.character(dt$SV_ANALYSIS_GROUP)
group_text <- ifelse(
  group_text %in% names(group_labels),
  group_labels[group_text],
  gsub("_", " ", group_text)
)
group <- factor(group_text, levels = unique(group_text))

row_ha <- rowAnnotation(
  `Gene score` = anno_barplot(
    dt$SCORE_NUM, border = FALSE, gp = gpar(fill = "grey45"),
    width = unit(18, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

col_fun <- colorRamp2(c(0, 1), c("#f2f2f2", "#303030"))
ht <- Heatmap(
  evidence,
  name = "Evidence",
  col = col_fun,
  na_col = "#d9d9d9",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = group,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 40,
  row_title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  row_title_rot = 0,
  row_gap = unit(3.5, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  left_annotation = row_ha,
  row_names_max_width = unit(62, "mm"),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("No", "Yes"),
    title = "Evidence\n(grey = unavailable)"
  ),
  column_title = "SV candidate evidence",
  column_title_gp = gpar(fontface = "bold", fontsize = 13)
)

height <- max(7, 0.30 * nrow(evidence) + 2.8)
pdf(pdf_file, width = 14.5, height = height, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(png_file, width = 14.5, height = height, units = "in", res = 400)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

svglite(svg_file, width = 14.5, height = height)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
