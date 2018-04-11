plot_cca_scores <- function(
  cca_obj, score = "xscores", cv = 1, n = 10
) {
  dat_xscore <- data.frame(
    y = extreme_n(cca_obj$scores[[score]][,cv], n = n)
  )
  dat_xscore$x <- rownames(dat_xscore)
  ggplot(dat_xscore, aes(x = reorder(x, y), y = y, fill = y)) +
    geom_segment(aes(xend = x, yend = 0), size = 0.2) +
    geom_point(shape = 21, size = 4, stroke = 0.2) +
    scale_fill_gradientn(colors = solar_flare) +
    geom_hline(yintercept = 0, size = 0.2) +
    coord_flip() +
    theme_clean(base_size = 20) +
    theme(legend.position = "none") +
    # labs(x = NULL, y = sprintf("xscores[,%i]", i))
    # labs(x = NULL, y = bquote("CV"["Bulk"~.(cv)]))
    labs(x = NULL, y = bquote("CV"[.(cv)]))
}