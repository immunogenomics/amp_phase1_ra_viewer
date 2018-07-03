plot_bulk_single_cca <- function(
  dat_cca, x = 1, y = 2, marker = ""
) {
  
  fill_values <- quantile_breaks(dat_cca$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  
  make_subplot <- function(
    this_data = "bulk", this_fill = "cell_type"
  ) {
    this_shape <- 22
    this_size <- 3
    if (this_data == "cell") {
      this_shape <- 21
      this_size <- 2
    }
    p <- ggplot() +
      geom_point(
        data = subset(
          dat_cca[order(dat_cca$marker),],
          data == this_data
        ),
        mapping = aes_string(
          x = sprintf("V%s", x),
          y = sprintf("V%s", y),
          fill = this_fill
        ),
        stroke = 0.1, shape = this_shape, size = this_size
      ) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      labs(x = sprintf("CV%s", x), y = sprintf("CV%s", y)) +
      theme_clean(base_size = 20) +
      theme(legend.position = "bottom")
    
    if (this_fill == "marker") {
      this_name <- bquote("Log"[2]~"(TPM+1)  ")
      if (this_data == "cell") {
        this_name <- bquote("Log"[2]~"(CPM+1)  ")
      }
      # scale_fill_manual(values = meta_colors$type) +
      p <- p + scale_fill_gradientn(
        colours = colorRampPalette(fill_palette)(length(fill_values)),
        values  = fill_values,
        breaks  = scales::pretty_breaks(n = 4),
        name    = this_name
      ) +
      guides(
        fill  = guide_colorbar(barwidth = 8, barheight = 1),
        # fill = guide_legend(
        #   nrow = 4, override.aes = list(shape = 22, size = 4)
        # ),
        shape = guide_legend(
          override.aes = list(stroke = 0.5, size = 4)
        )
      )
    } else if (this_fill == "cluster") {
      p <- p +
        scale_fill_manual(values = meta_colors$cluster, name = NULL) +
        guides(
          fill = guide_legend(
            nrow = 6,
            override.aes = list(shape = this_shape, size = 4)
          )
        )
    } else if (this_fill == "cell_type") {
      p <- p +
        scale_fill_manual(values = meta_colors$type, name = NULL) +
        guides(
          fill = guide_legend(
            ncol = 2,
            override.aes = list(shape = this_shape, size = 4)
          )
        )
    }
    return(p)
  }
  
  p1 <- make_subplot(this_data = "bulk", this_fill = "marker") +
    labs(subtitle = "bulk")
  
  p2 <- make_subplot(this_data = "cell", this_fill = "marker") +
    labs(subtitle = "cell")
  
  p3 <- make_subplot(this_data = "bulk", this_fill = "cell_type")
  
  p4 <- make_subplot(this_data = "cell", this_fill = "cluster")
  
  p1 + p2 + p3 + p4 + plot_layout(ncol = 2) +
    plot_annotation(
      title = marker,
      theme = theme(
        plot.title = element_text(
          size = 20, face = "italic", hjust = 0.5
        ),
        plot.caption = element_text(size = 20)
      )
    )
  
}

#plot_bulk_single_cca <- function(dat_cca, x = 1, y = 2) {
#  
#  fill_values <- quantile_breaks(dat_cca$marker, n = 9)
#  fill_values <- fill_values / max(fill_values)
#  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
#  
#  ggplot() +
#    # geom_hline(yintercept = 0, color = "grey90") +
#    # geom_vline(xintercept = 0, color = "grey90") +
#    geom_point(
#      data = dat_cca[order(dat_cca$marker),],
#      mapping = aes_string(
#        x = sprintf("V%s", x),
#        y = sprintf("V%s", y),
#        fill = "marker",
#        # fill = "cell_type",
#        shape = "data",
#        size = "data"
#      ),
#      stroke = 0.1
#    ) +
#    scale_shape_manual(values = c(22, 21)) +
#    scale_size_manual(values = c(3, 2)) +
#    # scale_fill_manual(values = meta_colors$type) +
#    scale_fill_gradientn(
#      colours = colorRampPalette(fill_palette)(length(fill_values)),
#      values  = fill_values,
#      breaks  = scales::pretty_breaks(n = 4),
#      name    = bquote("Log"[2]~"(CPM+1)  ")
#    ) +
#    guides(
#      fill  = guide_colorbar(barwidth = 10, barheight = 1),
#      # fill = guide_legend(
#      #   nrow = 4, override.aes = list(shape = 22, size = 4)
#      # ),
#      shape = guide_legend(
#        override.aes = list(stroke = 0.5, size = 4)
#      )
#    ) +
#    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
#    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
#    labs(x = sprintf("CV%s", x), y = sprintf("CV%s", y)) +
#    # geom_circle(
#    #   mapping = aes(x0 = 0, y0 = 0, r = 1)
#    # ) +
#    # coord_equal() +
#    facet_wrap(~ data, scales = "free") +
#    theme_clean(base_size = 20) +
#    theme(legend.position = "bottom")
#  
#}
