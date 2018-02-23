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
plot_tsne <- function(dat, tsne_x = "T1", tsne_y = "T2", title = NULL) {
  n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 2.0
  fill_values <- quantile_breaks(dat$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  theme_tsne <- theme_bw(base_size = 22) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(size = 25,  face="bold")
  )
  p1 <- ggplot() +
    geom_point(
      data    = dat[order(dat$marker),],
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "marker"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.15
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
    guides(
      fill  = guide_colorbar(barwidth = 10, barheight = 1),
      alpha = "none"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_tsne
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  p2 <- ggplot() +
    geom_point(
      data    = dat[sample(nrow(dat)),],
      # mapping = aes(x = T1, y = T2, fill = cluster),
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "cluster"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.15
    ) +
    # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
    guides(fill = guide_legend(nrow = 4, override.aes = list(size = 4))) +
    labs(x = NULL, y = NULL) +
    ggtitle("Identified clusters") +
    theme_tsne
  bottom_text <- sprintf(
    "%s cells, %s (%s%%) nonzero cells",
    nrow(dat),
    n_nonzero,
    signif(100 * n_nonzero / nrow(dat), 3)
  )
  egg::ggarrange(
    bottom = grid::textGrob(
      label = bottom_text, gp = grid::gpar(fontsize = 20)
    ),
    plots = list(p1, p2), ncol = 2
  )
}

