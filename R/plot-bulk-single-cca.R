plot_bulk_single_cca <- function(dat_cca, x = 1, y = 2) {
  fill_values <- quantile_breaks(dat_cca$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  ggplot() +
    # geom_hline(yintercept = 0, color = "grey90") +
    # geom_vline(xintercept = 0, color = "grey90") +
    geom_point(
      data = dat_cca,
      mapping = aes_string(
        x = sprintf("V%s", x),
        y = sprintf("V%s", y),
        fill = "marker",
        # fill = "cell_type",
        shape = "data",
        size = "data"
      ),
      stroke = 0.1
    ) +
    scale_shape_manual(values = c(22, 21)) +
    scale_size_manual(values = c(3, 2)) +
    # scale_fill_manual(values = meta_colors$type) +
    scale_fill_gradientn(
      colours = colorRampPalette(fill_palette)(length(fill_values)),
      values  = fill_values,
      breaks  = scales::pretty_breaks(n = 4),
      name    = bquote("Log"[2]~"(CPM+1)  ")
    ) +
    guides(
      fill  = guide_colorbar(barwidth = 10, barheight = 1),
      # fill = guide_legend(
      #   nrow = 4, override.aes = list(shape = 22, size = 4)
      # ),
      shape = guide_legend(
        override.aes = list(stroke = 0.5, size = 4)
      )
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(x = sprintf("CV%s", x), y = sprintf("CV%s", y)) +
    # geom_circle(
    #   mapping = aes(x0 = 0, y0 = 0, r = 1)
    # ) +
    # coord_equal() +
    facet_wrap(~ data, scales = "free") +
    theme_clean(base_size = 20) +
    theme(legend.position = "bottom")
}
