plot_bulk_dots <- function(dat, marker = "") {
  ggplot() +
    # geom_boxplot() +
    # geom_violin() +
    # geom_jitter(height = 0, width = 0.25, color = "dimgrey", size = 0.3) +
    # geom_quasirandom(
    # geom_sina(
    geom_jitter(
      data    = dat,
      mapping = aes(
        x = cell_type, y = marker, fill = inflamed, group = inflamed
      ),
      shape = 21, stroke = 0.15, size = 3.5,
      position = position_jitterdodge(jitter.width = 0.6)
    ) +
    geom_vline(xintercept = c(1,2,3) + 0.5, size = 0.1) +
    scale_y_continuous(breaks = scales::extended_breaks(n = 4)) +
    scale_x_discrete(limits = rev(levels(dat$cluster))) +
    geom_vline(
      data = data.frame(x = c(7.5, 11.5, 15.5)),
      mapping = aes(xintercept = x),
      color = "grey70"
    ) +
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
