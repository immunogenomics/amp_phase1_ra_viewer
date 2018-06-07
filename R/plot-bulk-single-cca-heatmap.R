plot_bulk_single_cca_heatmap <- function(dat_cca) {
  melt_cca <- data.table::melt(
    data = dat_cca,
    id.vars = c("cell_type", "cluster", "data")
  ) %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(scaled_value = scale(value)) %>%
    dplyr::group_by(variable, cluster) %>%
    dplyr::summarise(mean = mean(scaled_value))
  melt_cca <- seriate_columns(melt_cca, method = "BEA_TSP")
  ggplot() +
    geom_tile(
      # data = melt_cca,
      data = x,
      mapping = aes(x = variable, y = cluster, fill = mean),
      size = 0.2, color = "black"
    ) +
    coord_equal() +
    # scale_fill_viridis_c(
    #   breaks = quantile(x$mean, probs = seq(0, 1, length.out = 5))
    # ) +
    # scale_fill_gradientn(
    #   colors = solar_flare,
    #   name = "quantile",
    #   breaks = quantile(x$mean, probs = seq(0, 1, length.out = 5))
    # ) +
    scale_fill_gradientn(
      colors = solar_flare
    ) +
    scale_x_discrete(
      expand = c(0, 0), labels = function(x) substr(x, 2, 3)
    ) +
    guides(
      fill = guide_colorbar(barwidth = 0.8, barheight = 20)
    ) +
    labs(x = "CV", y = NULL) +
    # scale_y_discrete(expand = c(-0.05, 0.05)) +
    theme_minimal(base_size = 20) +
    theme(panel.grid = element_blank())
}
