plot_bulk_dots <- function(dat, marker = "") {
  
  # bulk.summary <- aggregate(marker ~ cell_type + inflamed, mean, data=dat)
  dat_median <- dat %>% group_by(cell_type) %>% group_by(inflamed) %>% summarise(median = median(marker))

  ggplot(data    = dat,
         mapping = aes(
           x = cell_type, y = dat$marker, fill = inflamed, group = inflamed
         )) +
    # geom_boxplot() +
    # geom_violin() +
    # geom_jitter(height = 0, width = 0.25, color = "dimgrey", size = 0.3) +
    # geom_quasirandom(
    # geom_sina(
    geom_jitter(
      shape = 21, stroke = 0.15, size = 3.5,
      position = position_jitterdodge(jitter.width = 0.2, seed = 1, dodge.width = 0.8)
    ) +
    # geom_segment(data = bulk.summary, aes(x=cell_type,xend=max(pos),y=0,yend=0),colour="green") + 
    stat_summary(
      fun.data=mean_sdl, fun.args = list(mult=1),
      geom="pointrange", 
      size = 0.6,
      shape = 124,
      fatten = 8,
      position = position_jitterdodge(jitter.width = 0.2, seed = 1, dodge.width = 0.8),
      color="black",
      show.legend = FALSE
      ) +
    # stat_summary(
    #   fun.y = median, fun.ymin = median, fun.ymax = median,
    #   geom = "crossbar",
    #   position = position_jitterdodge(jitter.width = 0.2, seed = 1, dodge.width = 0.8),
    #   width = 0.8, size = 0.5,
    #   colour = "black",
    #   show.legend = FALSE
    #   # colour = rep(c("#FF7F00", "#FFD8B2", "#6A3D9A"), 4)
    # ) +
    # geom_boxplot(position = position_dodge(.75)) + 
    # geom_hline(data = bulk.summary, aes(yintercept = marker, width = 0.3), 
    #            col = rep(c("#6A3D9A", "#FFD8B2", "#FF7F00"), 4)) +
    # geom_crossbar(data=bulk.summary,aes(x = c(0.5, 0, -0.5) + cell_type, 
    #                                     y = marker,
    #                                     ymax = marker,
    #                                     ymin = marker),
    #                size = 0.3, col = rep(c("#6A3D9A", "#FFD8B2", "#FF7F00"), 4), width = 0.3) +
    geom_vline(xintercept = c(1,2,3) + 0.5, size = 0.1) +
    scale_y_continuous(breaks = scales::extended_breaks(n = 4)) +
    scale_x_discrete(limits = rev(levels(dat$cell_type))) +
    # geom_vline(
    #   data = data.frame(x = c(7.5, 11.5, 15.5)),
    #   mapping = aes(xintercept = x),
    #   color = "grey70"
    # ) +
    coord_flip() +
    labs(
      x     = NULL,
      y     = bquote("Log"[2]~"(TPM+1)"),
      title = marker
      # subtitle = tsne_subtitle
    ) +
    scale_fill_manual(
      values = meta_colors$inflamed, name = NULL
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(size = 4),
        reverse = TRUE
      )
    ) +
    theme_bw(base_size = 26) + theme(
      # legend.position = "none",
      legend.key.size = unit(1.5, "lines"),
      axis.ticks      = element_line(size = 0.5),
      # axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      # axis.text       = element_blank(),
      # axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title      = element_text(size = 26, face = "italic")
    )
}
