plot_box <- function(dat, title = "") {
  
  # Put the clusters in the desired order.
  clusters <- c(
    "SC-M1",  "SC-M2", "SC-M3", "SC-M4",
    "SC-F1", "SC-F2", "SC-F3", "SC-F4",
    "SC-B1", "SC-B2", "SC-B3", "SC-B4",
    "SC-T1", "SC-T2", "SC-T3", "SC-T4", "SC-T5", "SC-T6"
  )
  
  # ix_zero <- which(dat$marker == 0)
  # dat$marker[ix_zero] <- dat$marker[ix_zero] + runif(n = length(ix_zero), min = -0.5, max = 0)
  
  dat$cluster <- factor(dat$cluster, levels = clusters)
  p1 <- ggplot() +
    # geom_boxplot() +
    # geom_violin() +
    # geom_jitter(height = 0, width = 0.25, color = "dimgrey", size = 0.3) +
    geom_quasirandom(
      data    = subset(dat, marker > 0),
      # data    = dat,
      mapping = aes(x = cluster, y = marker, fill = cluster),
      shape = 21, stroke = 0.15, size = 2.5
    ) +
    scale_y_continuous(breaks = scales::extended_breaks(n = 4)) +
    scale_x_discrete(limits = rev(levels(dat$cluster))) +
    geom_vline(
      data = data.frame(x = c(6.5, 10.5, 14.5)),
      mapping = aes(xintercept = x),
      color = "grey70"
    ) +
    coord_flip() +
    labs(
      x     = NULL,
      y     = bquote("Log"[2]~"(CPM+1)")
      # title = title
      # subtitle = tsne_subtitle
    ) +
    scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    theme_bw(base_size = 26) + theme(
      legend.position = "none",
      axis.ticks      = element_line(size = 0.5),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      # axis.text       = element_blank(),
      # axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title = element_text(size = 30,  face="bold")
    )
  # p1
  
  dat_percent <- dat %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(percent = sum(marker > 0) / length(marker) * 100)
  
  # Add small columns on the left showing percent of cells with 0 expression.
  # data_cols <- data.frame(
  #   x = 19:1, y = -1, cluster = NA
  # )
  # p1 <- p1 + geom_col(
  #     data = dat_percent,
  #     mapping = aes(x = cluster, y = -1 * (100 - percent) / 100, fill = cluster),
  #     size = 0.15, colour = "black", width = 0.5
  #   ) + geom_col(
  #     data = data_cols,
  #     mapping = aes(x = x, y = y, fill = cluster),
  #     size = 0.15, colour = "black", width = 0.5
  #   )
  
  p2 <- ggplot() +
    geom_col(
      data = dat_percent,
      mapping = aes(x = cluster, y = percent, fill = cluster),
      color = "black", size = 0.15, width = 0.8
      # fill = "grey60"
    ) +
    geom_vline(
      data = data.frame(x = c(6.5, 10.5, 14.5)),
      mapping = aes(xintercept = x),
      color = "grey70"
    ) +
    scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    scale_x_discrete(limits = rev(levels(dat$cluster))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    labs(x = NULL, y = "% non-zero", title = title) +
    coord_flip() +
    theme_bw(base_size = 26) + theme(
      legend.position = "none",
      axis.ticks      = element_line(size = 0.5),
      # axis.text.y       = element_blank(),
      # axis.ticks.y      = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title      = element_text(size = 25, face = "italic")
    )
  p2 + p1 + plot_layout(widths = c(0.4, 0.6))
  
  # proportion <- rep(0, length(table(dat$cluster)))
  # for (i in 1:length(table(dat$cluster))){
  #   proportion[i] <- sum(dat$cluster == i & dat$marker > 0)/ (table(dat$cluster)[i])
  #   # print(sum(dat$cluster == i & dat$marker > 0))
  #   # print(table(dat$cluster)[i])
  # }
  # dat_pro <- data.frame(
  #   cluster = as.character(seq(1, length(table(dat$cluster)))),
  #   nonzero = proportion
  # )
  # p2 <- ggplot(
  #   data=dat_pro, 
  #   aes(x=cluster, y= nonzero, fill = cluster)
  #   # aes(x=cluster, y= percent(nonzero), fill = cluster)
  #   ) +
  #   geom_bar(stat="identity", position = "stack") +
  #   labs(
  #     x = NULL,
  #     y = "Percent nonzero",
  #     title = title
  #   ) +
  #   scale_y_continuous(labels = percent) +
  #   # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
  #   scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
  #   theme_box
  # p3 <- ggplot()
  # bottom_text <- sprintf(
  #   "%s is expressing the highest percent of nonzero cells in cluster %s.",
  #   title,
  #   as.integer(dat_pro$cluster[which(dat_pro$nonzero == max(dat_pro$nonzero))])
  # )
}

