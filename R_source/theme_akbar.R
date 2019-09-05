.theme.akbar<-function(base_size = 14, base_family = ""){
  #mat's ggplot2 theme http://doi.org/10.5281/zenodo.3362026
  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
     ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5),
                    strip.background = ggplot2::element_rect(fill = grDevices::rgb(0,20, 100, 30, maxColorValue = 100), linetype = 0),
                    panel.background = ggplot2::element_rect(fill = grDevices::rgb(0, 0, 100, 3, maxColorValue = 100), linetype = 0, color = NA),
                    line = ggplot2::element_blank(),
                    text = ggplot2::element_text(),
                    plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                    panel.border = ggplot2::element_rect(colour = NA),
                    axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
                    axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
                    axis.title.x = ggplot2::element_text(vjust = -0.2), axis.text = ggplot2::element_text(),
                    axis.line = ggplot2::element_line(colour = "black"), axis.ticks = ggplot2::element_line(),
                    panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                    legend.key = ggplot2::element_blank(), legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.key.size = grid::unit(0.2, "cm"),
                    legend.title = ggplot2::element_text(face = "italic", vjust = 1),
                    legend.box.background = ggplot2::element_blank(),
                    legend.background = ggplot2::element_rect(fill = "transparent",color = NA),
                    plot.margin = grid::unit(c(10, 5, 5, 5),"mm"),
                    strip.text = ggplot2::element_text(face = "bold")))
}
