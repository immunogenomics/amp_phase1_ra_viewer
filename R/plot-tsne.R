#' Create 2 tSNE plots side by side.
#' The left plot is colored by marker.
#' The right plot is colored by cluster.
#' @param dat A dataframe with columns T1, T2, marker, cluster
plot_tsne <- function(dat, tsne_x = "T1", tsne_y = "T2", title = NULL) {
  n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 2.0
  if (nrow(dat) < 5000) {
    point_size <- 2.5
  }
  fill_values <- quantile_breaks(dat$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  theme_tsne_1 <- theme_bw(base_size = 25) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 30),
    legend.text    = element_text(size = 18)
  )
  theme_tsne_2 <- theme_bw(base_size = 25) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 30),
    legend.text    = element_text(size = 18)
  )
  p1 <- ggplot() +
    geom_point(
      data    = dat[order(dat$marker),],
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
    guides(
      fill  = guide_colorbar(
        ticks.linewidth = 1,
        ticks.colour = "black",
        barwidth = 7,
        barheight = 1.7
      ),
      alpha = "none"
    ) +
    labs(x = NULL, y = NULL, title = substitute(italic(x), list(x = title))) +
    theme_tsne_1
  
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  p2 <- ggplot() +
    geom_point(
      data    = dat[sample(nrow(dat)),],
      # mapping = aes(x = T1, y = T2, fill = cluster),
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "cluster"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.1
    ) +
    # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    guides(
      fill = guide_legend(
        nrow = 6,
        override.aes = list(size = 6)
      )
    ) +
    labs(x = NULL, y = NULL, title = "Clusters") +
    # ggtitle("Identified clusters") +
    theme_tsne_2
  
  bottom_text <- sprintf(
    "%s cells: %s (%s%%) non-zero",
    comma(nrow(dat)),
    comma(n_nonzero),
    signif(100 * n_nonzero / nrow(dat), 2)
  )
  
  p1 + p2 + plot_annotation(
    caption = bottom_text,
    theme = theme(
      plot.caption = element_text(size = 25, hjust = 0.5)
    )
  )
}

# marker <- "THY1"
# gene_ix <- which(gene_symbols == marker)
# meta$marker <- as.numeric(log2cpm[gene_ix,])
# cell_ix <- seq(nrow(meta))
# tsne_x <- "T1_all"
# tsne_y <- "T2_all"
# plot_tsne(meta[cell_ix,], tsne_x, tsne_y, title = marker)
