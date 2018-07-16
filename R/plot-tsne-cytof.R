#' Get the quantile breaks in a numeric vector.
#' @param x A numeric vector.
#' @param n The number of breaks.
#' @return A vector with unique breaks.
quantile_breaks <- function(x, n = 10) {
  breaks <- quantile(x, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#' Create 2 tSNE plots side by side.
#' The left plot is colored by marker.
#' The right plot is colored by cluster.
#' @param dat A dataframe with columns T1, T2, marker, cluster
plot_tsne_cytof <- function(dat, tsne_x = "T1", tsne_y = "T2", title = NULL) {
  # n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 1.0
  if (nrow(dat) < 5000) {
    point_size <- 1.5
  }
  fill_values <- quantile_breaks(dat$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  theme_tsne <- theme_bw(base_size = 15) + theme(
    legend.position = "right",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 22),
    legend.text    = element_text(size = 15)
  )
  p1 <- ggplot() +
    geom_point(
      # data    = dat[order(dat$marker),],
      data = dat,
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "marker"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.1
    ) +
    scale_fill_gradientn(
      # Linear scale
      # colours = fill_palette,
      # Quantile scale
      colours = colorRampPalette(fill_palette)(length(fill_values)),
      values  = fill_values,
      breaks  = scales::pretty_breaks(n = 4),
      name    = bquote("Log"[2]~"(CPM+1)  ")
    ) +
    labs(x = NULL, y = NULL, title = substitute(italic(x), list(x = title))) +
    guides(
    fill  = guide_colorbar(barwidth = 7, barheight = 0.7),
    alpha = "none"
    ) +
    theme_tsne
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  p2 <- ggplot() +
    geom_point(
      # data    = dat[sample(nrow(dat)),],
      data = dat,
      # mapping = aes(x = T1, y = T2, fill = cluster),
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "cluster"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.1
    ) +
    # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = meta_colors$cytof_cluster, name = "Cluster") +
    labs(x = NULL, y = NULL, title = "Clusters") +
    guides(fill = guide_legend(nrow = 10)) + # , override.aes = list(size = 2)
    theme_tsne
  # bottom_text <- sprintf(
  #   "%s cells, %s (%s%%) nonzero cells",
  #   nrow(dat),
  #   n_nonzero,
  #   signif(100 * n_nonzero / nrow(dat), 3)
  # )
  p1 + p2 + plot_annotation(
    # caption = bottom_text,
    theme = theme(
      plot.caption = element_text(size = 18, hjust = 0.5)
    )
  )
}

