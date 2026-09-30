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

dt <- fread(input_file, na.strings = c("", ".", "NA", "N/A"))
if (nrow(dt) == 0) stop("No SV-gene candidates available for heatmap.")

dt[, CALLER_COUNT_NUM := suppressWarnings(as.numeric(CALLER_COUNT))]
dt[is.na(CALLER_COUNT_NUM), CALLER_COUNT_NUM := 0]
dt[, SCORE_NUM := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
dt[is.na(SCORE_NUM), SCORE_NUM := 0]
setorder(dt, -SCORE_NUM, -CALLER_COUNT_NUM)
dt <- unique(dt, by = c("SV_ID", "GENE"))
dt <- head(dt, top_n)

flag <- function(x, positive) {
  ifelse(is.na(x), NA_real_, ifelse(x %in% positive, 1, 0))
}

effect_direct <- !grepl("INSIDE_INVERSION|UNRESOLVED", dt$SV_GENE_EFFECT)
evidence <- cbind(
  flag(dt$CALL_SUPPORT, "MULTI_CALLER"),
  flag(dt$POPULATION_STATUS, "RARE"),
  flag(dt$PANEL_STATUS, "PANEL_GENE"),
  flag(dt$MITOCARTA, c("NUCLEAR_MITOCHONDRIAL_GENE", "MTDNA_ENCODED_GENE")),
  ifelse(effect_direct, 1, 0),
  flag(dt$LONGPHASE_PHASED, "YES"),
  flag(dt$STRAGLR, "YES"),
  flag(dt$TLDR, "YES"),
  flag(dt$METHYLATION_CONTEXT, "EVALUATED")
)
colnames(evidence) <- c(
  "Multi-caller", "Rare in needLR", "ON panel gene", "MitoCarta gene",
  "Direct/near gene effect", "LongPhase phased", "Straglr match",
  "TLDR match", "Methylation evaluated"
)

rownames(evidence) <- paste0(
  dt$GENE, " | ", dt$SVTYPE, " | ", dt$CHROM, ":",
  format(as.numeric(dt$START), scientific = FALSE, trim = TRUE)
)
group <- factor(dt$SV_ANALYSIS_GROUP)

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
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 35,
  rect_gp = gpar(col = "white", lwd = 0.8),
  left_annotation = row_ha,
  heatmap_legend_param = list(at = c(0, 1), labels = c("No", "Yes")),
  column_title = "SV candidate evidence",
  column_title_gp = gpar(fontface = "bold", fontsize = 13)
)

height <- max(7, 0.30 * nrow(evidence) + 2.8)
pdf(pdf_file, width = 12.5, height = height, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(png_file, width = 12.5, height = height, units = "in", res = 400)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

svglite(svg_file, width = 12.5, height = height)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
