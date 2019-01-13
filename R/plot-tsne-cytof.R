#' Get the quantile breaks in a numeric vector.
#' @param x A numeric vector.
#' @param n The number of breaks.
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
      name    = "Normalized\nIntensity"
    ) +
    labs(x = NULL, y = NULL, title = substitute(x, list(x = title))) +
    guides(
      fill  = guide_colorbar(
        ticks.linewidth = 1,
        ticks.colour = "black",
        barwidth = 8,
        barheight = 1.7
        # frame.colour = "black"
      ),
      alpha = "none"
    ) +
    theme_tsne_1
  
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  p2 <- ggplot() +
    geom_point(
      data    = dat[sample(nrow(dat)),],
      mapping = aes_string(x = tsne_x, y = tsne_y, fill = "cluster"),
      size    = point_size,
      shape   = 21,
      stroke  = 0.1
    ) +
    scale_fill_manual(values = meta_colors$cytof_cluster, name = "Cluster") +
    labs(x = NULL, y = NULL, title = "Clusters") +
    guides(
      fill = guide_legend(nrow = 10, override.aes = list(size = 6))
    ) +
    theme_tsne_2
  
  n_nonzero <- sum(dat$marker > 0)
  bottom_text <- sprintf(
    "%s cells: %s (%s%%) non-zero",
    comma(nrow(dat)),
    comma(n_nonzero),
    signif(100 * n_nonzero / nrow(dat), 2)
  )
  
  p1 + p2 + plot_annotation(
    caption = bottom_text,
    theme = theme(
      plot.caption = element_text(size = 20, hjust = 0.5)
    )
  )
}

# marker <- "CD90"
# cytof_all$marker <- as.numeric(cytof_all[, which(colnames(cytof_all) == marker)])
# cell_ix <- which(cytof_all$cell_type == "Fibroblast")
# plot_tsne_cytof(cytof_all[cell_ix,], "SNE1", "SNE2", title = marker)

# cytof_summarize$pct_nonzero <- as.numeric(substr(cytof_summarize$pct_nonzero, 1, 5))
# 
# cytof_dat <- cytof_all %>%
#   dplyr::select(-sampleID, -site, -SNE1, -SNE2, -status.tissue, -type, -status) %>%
#   data.table::melt(id.vars = c("cell_type", "cluster")) %>%
#   dplyr::group_by(cell_type, cluster, variable) %>%
#   dplyr::summarize(pct_nonzero = sqrt(mean(value, na.rm = TRUE)))
#   # dplyr::summarize(pct_nonzero = sum(value > 0) / length(value))
# cytof_dat$cluster <- as.character(cytof_dat$cluster)
# cytof_dat$variable <- as.character(cytof_dat$variable)
# 
# cytof_mat <- data.table::dcast(
#   data = cytof_dat,
#   formula = cell_type + cluster ~ variable,
#   value.var = "pct_nonzero"
# )
# 
# cytof_ord <- seriation::seriate(
#   x = as.matrix(cytof_mat[,!colnames(cytof_mat) %in% c("cell_type", "cluster")]),
#   method = "PCA"
# )
# cytof_mat <- cytof_mat[cytof_ord[[1]], c(1, 2, cytof_ord[[2]])]
# cytof_mat[1:5,1:5]                   
# 
# # cytof_dat$cluster <- factor(cytof_dat$cluster, cytof_mat$cluster[cytof_ord[[1]]])
# # cytof_dat$variable <- factor(
# #   x = cytof_dat$variable, 
# #   levels = colnames(cytof_mat[,!colnames(cytof_mat) %in% c("cell_type", "cluster")])[cytof_ord[[2]]]
# # )
# 
# ggplot(cytof_dat) +
#   geom_tile(aes(x = variable, y = cluster, fill = pct_nonzero)) +
#   scale_fill_gradient2(low = "grey90", high = "red", na.value = "white") +
#   facet_grid(cell_type ~ ., scales = "free_y", space = "free_y")
