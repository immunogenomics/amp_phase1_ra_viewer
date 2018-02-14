plot_box <- function(dat, title = "") {
  
  # Put the clusters in the desired order.
  clusters <- c(
    "CB1", "CB2", "CB3", "CB4",
    "CT1", "CT2", "CT3", "CT4", "CT5", "CT6", "CT7",
    "CM1",  "CM2", "CM3", "CM4",
    "CF1", "CF2", "CF3", "CF4"
  )
  
  dat$cluster <- factor(dat$cluster, levels = clusters)
  p1 <- ggplot(
    data    = subset(dat, marker > 0),
    mapping = aes(x = cluster, y = marker, fill = cluster)) +
    # geom_boxplot() +
    # geom_violin() +
    # geom_jitter(height = 0, width = 0.25, color = "dimgrey", size = 0.3) + 
    geom_quasirandom(shape = 21, stroke = 0.1, size = 2) +
    scale_x_discrete(limits = rev(levels(dat$cluster))) +
    coord_flip() +
    labs(
      x     = NULL,
      y     = bquote("Log"[2]~"(CPM+1)"),
      title = title
      # subtitle = tsne_subtitle
    ) +
    scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
    theme_bw(base_size = 22) + theme(
      legend.position = "none",
      # axis.text       = element_blank(),
      # axis.ticks      = element_blank(),
      # panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title = element_text(size = 25,  face="bold")
    )
  p1
  
  dat_percent <- dat %>%
    group_by(cluster) %>%
    summarise(percent = sum(marker > 0) / length(marker) * 100)
  
  p2 <- ggplot() +
    geom_col(
      data = dat_percent,
      mapping = aes(x = cluster, y = percent, fill = cluster),
      color = "grey50", size = 0.1
      # fill = "grey60"
    ) +
    scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
    scale_x_discrete(limits = rev(levels(dat$cluster))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    labs(x = NULL, y = "% non-zero") +
    coord_flip() +
    theme_bw(base_size = 22) + theme(
      legend.position = "none",
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      # panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title = element_text(size = 25,  face="bold")
    )
  
  egg::ggarrange(p1, p2, ncol = 2, widths = c(0.8, 0.2))
  
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
  #   scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  #   theme_box
  # p3 <- ggplot()
  # bottom_text <- sprintf(
  #   "%s is expressing the highest percent of nonzero cells in cluster %s.",
  #   title,
  #   as.integer(dat_pro$cluster[which(dat_pro$nonzero == max(dat_pro$nonzero))])
  # )
  # egg::ggarrange(
  #   # bottom = textGrob(
  #   #   label = bottom_text, 
  #   #   gp = gpar(fontsize = 20, fontface="bold")
  #   # ),
  #   plots = list(p1, p2), ncol = 2, widths = 2.7:1, height = 2.2:1,
  #   hr()
  # )
}

