plot_bulk_dots <- function(dat, marker = "") {
  ps <- sapply(unique(dat$cell_type), function(i) {
    kruskal.test(formula = marker ~ inflamed, data = subset(dat, cell_type == i))$p.value
  })
  ps <- sprintf("p = %s", format.pval(ps, 1))
  ps[2] <- sprintf("Kruskal-Wallis, %s", ps[2])
  dat_p <- data.frame(
    cell_type = factor(
      x = unique(as.character(dat$cell_type)),
      levels = levels(dat$cell_type)
    ),
    pvalue = ps
  )
  dat$inflamed <- factor(
    x = as.character(dat$inflamed),
    levels = c("blank1", levels(dat$inflamed))
  )
  # dat <- b_meta
  ggplot() +
  geom_boxplot(
    data    = dat,
    mapping = aes(
      x = inflamed, y = marker, fill = inflamed, group = inflamed
    ),
    color = "grey20",
    alpha = 0.4, size = 0.4, width = 0.8, outlier.alpha = 0
  ) +
  geom_point(
    data    = dat,
    mapping = aes(
      x = inflamed, y = marker, fill = inflamed, group = inflamed
    ),
    shape = 21, size = 2.5,
    # position = position_jitterdodge(jitter.width = 0.9)
    position = position_quasirandom()
  ) +
  facet_grid(cell_type ~ ., switch = "both") +
  scale_x_discrete(
    limits = c("OA", "leukocyte-poor RA", "leukocyte-rich RA", "blank1")
  ) +
  scale_y_continuous(
    breaks = scales::extended_breaks(n = 4)
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
  geom_text(
    data = dat_p, size = 6,
    mapping = aes(x = 4, y = 0, hjust = 0, label = pvalue)
  ) +
  guides(
    fill = guide_legend(keyheight = unit(3, "lines"),
      override.aes = list(size = 0.5, alpha = 0.9),
      reverse = TRUE
    )
  ) +
  theme_bw(base_size = 26) + theme(panel.spacing.y = unit(0.5, "lines"),
    strip.text.y = element_text(angle = 180, hjust = 1),
    strip.background = element_blank(),
    legend.box.spacing = unit(0.1, "lines"),
    # legend.position = "none",
    legend.key.size = unit(1.5, "lines"),
    axis.ticks      = element_line(size = 0.5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    # axis.text       = element_blank(),
    # axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 26, face = "italic")
  )
}
# marker <- "THY1"
# b_meta$marker <- as.numeric(b_log2tpm[marker,])
# plot_bulk_dots(b_meta, marker)

