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
score_col <- if ("FINAL_GENE_RELEVANCE_DISPLAY_SCORE" %in% names(dt)) {
  "FINAL_GENE_RELEVANCE_DISPLAY_SCORE"
} else {
  "GENE_RELEVANCE_SCORE"
}
dt[, SCORE_NUM := suppressWarnings(as.numeric(get(score_col)))]
dt[is.na(SCORE_NUM), SCORE_NUM := 0]
dt[, RARE_FLAG := as.integer(
  POPULATION_CLASS %in% c(
    "LOW_FREQUENCY_BY_NEEDLR",
    "LOW_FREQUENCY_BY_EXACT_GNOMAD",
    "LOW_FREQUENCY_SUPPORTED_BY_NEEDLR_AND_GNOMAD"
  )
)]
dt[, PANEL_FLAG := as.integer(PANEL_STATUS == "PANEL_GENE")]

# Global overview is computed before top-gene selection so every gene in the
# candidate table contributes to the population-by-effect summary.
needlr_labels <- c(
  "NOT_OBSERVED_IN_NEEDLR_CONTROLS" = "Not observed",
  "VERY_RARE_NEEDLR_LE_0.001" = "Very rare <=0.1%",
  "RARE_NEEDLR_LE_0.01" = "Rare <=1%",
  "COMMON_NEEDLR_GT_0.01" = "Common >1%",
  "VERY_COMMON_NEEDLR_GE_0.05" = "Very common >=5%",
  "NO_NEEDLR_MATCH_AF_UNKNOWN" = "No match / AF unknown",
  "NEEDLR_NOT_EVALUABLE_GE_10MB" = "Not evaluable >=10 Mb",
  "NEEDLR_NOT_EVALUABLE_BREAKEND" = "Breakend not evaluable",
  "NEEDLR_FREQUENCY_UNKNOWN" = "Frequency unknown"
)

gnomad_labels <- c(
  "EXACT_GNOMAD_AF_ZERO" = "Exact match AF=0",
  "EXACT_GNOMAD_VERY_RARE_LE_0.001" = "Exact match <=0.1%",
  "EXACT_GNOMAD_RARE_LE_0.01" = "Exact match <=1%",
  "EXACT_GNOMAD_COMMON_GT_0.01" = "Exact match >1%",
  "EXACT_GNOMAD_VERY_COMMON_GE_0.05" = "Exact match >=5%",
  "EXACT_GNOMAD_MATCH_AF_UNAVAILABLE" = "Exact match / AF unavailable",
  "NO_EXACT_GNOMAD_SV_MATCH" = "No exact match",
  "GNOMAD_SV_RESOURCE_NOT_CONFIGURED" = "Resource not configured"
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

dt[, NEEDLR_LABEL := fifelse(
  NEEDLR_FREQUENCY_CLASS %in% names(needlr_labels),
  needlr_labels[NEEDLR_FREQUENCY_CLASS],
  gsub("_", " ", NEEDLR_FREQUENCY_CLASS)
)]
dt[, GNOMAD_LABEL := fifelse(
  GNOMAD_EXACT_AF_CLASS %in% names(gnomad_labels),
  gnomad_labels[GNOMAD_EXACT_AF_CLASS],
  gsub("_", " ", GNOMAD_EXACT_AF_CLASS)
)]
dt[, EFFECT_LABEL := fifelse(
  SV_EFFECT_GROUP %in% names(effect_labels),
  effect_labels[SV_EFFECT_GROUP],
  gsub("_", " ", SV_EFFECT_GROUP)
)]

global <- unique(dt[, .(GENE, NEEDLR_LABEL, EFFECT_LABEL)])
global <- global[, .(UNIQUE_GENES = uniqueN(GENE)), by = .(NEEDLR_LABEL, EFFECT_LABEL)]

p0 <- ggplot(global, aes(EFFECT_LABEL, NEEDLR_LABEL, fill = UNIQUE_GENES)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = ifelse(UNIQUE_GENES > 0, UNIQUE_GENES, "")), size = 3) +
  scale_fill_gradient(low = "grey95", high = "grey20") +
  labs(
    title = "All SV-overlapping genes",
    subtitle = "Each cell counts unique genes with at least one SV in the corresponding population-frequency and SV-effect classes.",
    x = "SV-gene relationship",
    y = "Population-frequency context",
    fill = "Unique genes"
  ) +
  theme_thesis(9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    plot.subtitle = element_text(size = 8.5, margin = margin(b = 6))
  )

rank_col <- if ("FINAL_EVENT_RANK_WITHIN_PANEL_STATUS" %in% names(dt)) {
  "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS"
} else if ("EVENT_RANK_WITHIN_PANEL_STATUS" %in% names(dt)) {
  "EVENT_RANK_WITHIN_PANEL_STATUS"
} else {
  NA_character_
}

dt[, PANEL_GROUP := fifelse(PANEL_STATUS == "PANEL_GENE", "Panel", "Non-panel")]
gene_rank <- dt[, .(
  RARE_COUNT = sum(RARE_FLAG),
  SCORE = max(SCORE_NUM),
  SV_COUNT = uniqueN(SV_ID),
  BEST_GROUP_RANK = if (!is.na(rank_col)) {
    suppressWarnings(min(as.numeric(get(rank_col)), na.rm = TRUE))
  } else {
    NA_real_
  }
), by = .(PANEL_GROUP, GENE)]
gene_rank[!is.finite(BEST_GROUP_RANK), BEST_GROUP_RANK := NA_real_]
gene_rank[, RANK_MISSING := as.integer(is.na(BEST_GROUP_RANK))]
setorder(
  gene_rank,
  PANEL_GROUP,
  RANK_MISSING,
  BEST_GROUP_RANK,
  -SCORE,
  -RARE_COUNT,
  -SV_COUNT,
  GENE
)
selected_table <- gene_rank[, head(.SD, top_n), by = PANEL_GROUP]
selected <- selected_table$GENE
plot_dt <- dt[GENE %in% selected]
ordered_genes <- c(
  selected_table[PANEL_GROUP == "Panel", GENE],
  selected_table[PANEL_GROUP == "Non-panel", GENE]
)
plot_dt[, GENE := factor(GENE, levels = rev(unique(ordered_genes)))]


needlr <- plot_dt[, .(N = uniqueN(SV_ID)), by = .(GENE, NEEDLR_LABEL)]
gnomad <- plot_dt[, .(N = uniqueN(SV_ID)), by = .(GENE, GNOMAD_LABEL)]
effect <- plot_dt[, .(N = uniqueN(SV_ID)), by = .(GENE, EFFECT_LABEL)]

p1 <- ggplot(needlr, aes(NEEDLR_LABEL, GENE, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(N > 0, N, "")), size = 2.5) +
  scale_fill_gradient(low = "grey95", high = "grey20") +
  labs(
    title = "needLR frequency context",
    x = NULL, y = NULL, fill = "SV count"
  ) +
  theme_thesis(9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right"
  )

p2 <- ggplot(gnomad, aes(GNOMAD_LABEL, GENE, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(N > 0, N, "")), size = 2.5) +
  scale_fill_gradient(low = "grey95", high = "grey20") +
  labs(
    title = "Exact gnomAD-SV v4.1 frequency context",
    x = NULL, y = NULL, fill = "SV count"
  ) +
  theme_thesis(9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

p3 <- ggplot(effect, aes(EFFECT_LABEL, GENE, fill = N)) +
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

combined <- p0 / (p1 + p2 + p3) +
  plot_layout(heights = c(1.0, 1.7)) +
  plot_annotation(
    title = "Population context and structural-variant effects by gene",
    subtitle = paste(
      "needLR and gnomAD-SV frequencies are displayed separately. gnomAD AF is reported only for the conservative exact coordinate/type match;",
      "AnnotSV benign-region overlap remains a separate context field and is not substituted for gnomAD exact-site AF."
    )
  )

height <- max(11, 0.28 * length(selected) + 6.0)
ggsave(pdf_file, combined, width = 19, height = height, device = cairo_pdf)
ggsave(png_file, combined, width = 19, height = height, dpi = 400)
ggsave(svg_file, combined, width = 19, height = height, device = svglite)
