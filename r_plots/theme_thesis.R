library(ggplot2)

theme_thesis <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 12, 8, 8)
    )
}
