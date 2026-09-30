library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)

phenotype_file <- snakemake@input[["phenotypes"]]
genes_file <- snakemake@input[["genes"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]
svg_file <- snakemake@output[["svg"]]
top_genes <- as.integer(snakemake@params[["top_genes"]])
top_hpo <- as.integer(snakemake@params[["top_hpo"]])

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

pheno <- fread(phenotype_file, na.strings = c("", ".", "NA"))
genes <- fread(genes_file, na.strings = c("", ".", "NA"))
genes[, SCORE := suppressWarnings(as.numeric(GENE_RELEVANCE_SCORE))]
genes[is.na(SCORE), SCORE := 0]
selected_genes <- head(genes[order(-SCORE)]$GENE, top_genes)

pheno <- pheno[gene_symbol %in% selected_genes & grepl("^HP:", hpo_id)]
if (nrow(pheno) == 0) {
  pdf(pdf_file, width = 10, height = 4.5, useDingbats = FALSE)
  grid.newpage()
  grid.text("No HPO annotations available for selected genes", gp = gpar(fontsize = 14))
  dev.off()
  png(png_file, width = 10, height = 4.5, units = "in", res = 400)
  grid.newpage()
  grid.text("No HPO annotations available for selected genes", gp = gpar(fontsize = 14))
  dev.off()
  svglite(svg_file, width = 10, height = 4.5)
  grid.newpage()
  grid.text("No HPO annotations available for selected genes", gp = gpar(fontsize = 14))
  dev.off()
  quit(save = "no", status = 0)
}

term_counts <- pheno[, .(N = uniqueN(gene_symbol)), by = .(hpo_id, hpo_label)]
selected_terms <- head(term_counts[order(-N, hpo_id)], top_hpo)
pheno <- pheno[hpo_id %in% selected_terms$hpo_id]

matrix_dt <- unique(pheno[, .(gene_symbol, hpo_id)])
matrix_dt[, value := 1L]
wide <- dcast(matrix_dt, gene_symbol ~ hpo_id, value.var = "value", fill = 0)
gene_order <- selected_genes[selected_genes %in% wide$gene_symbol]
wide <- wide[match(gene_order, gene_symbol)]
mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$gene_symbol

term_names <- selected_terms$hpo_label
names(term_names) <- selected_terms$hpo_id
colnames(mat) <- paste0(colnames(mat), "\n", term_names[colnames(mat)])

gene_info <- genes[match(rownames(mat), GENE)]
row_ha <- rowAnnotation(
  Panel = ifelse(gene_info$PANEL_STATUS == "PANEL_GENE", "Panel", "Non-panel"),
  Score = anno_barplot(
    gene_info$SCORE, border = FALSE, gp = gpar(fill = "grey45"),
    width = unit(18, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

ht <- Heatmap(
  mat,
  name = "HPO",
  col = c("0" = "#f5f5f5", "1" = "#303030"),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 7.5),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.8),
  left_annotation = row_ha,
  column_title = "Human HPO associations of prioritized genes",
  column_title_gp = gpar(fontface = "bold", fontsize = 13),
  heatmap_legend_param = list(at = c(0, 1), labels = c("No annotation", "Annotated"))
)

height <- max(7, 0.32 * nrow(mat) + 2.8)
pdf(pdf_file, width = 14, height = height, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(png_file, width = 14, height = height, units = "in", res = 400)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
svglite(svg_file, width = 14, height = height)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
