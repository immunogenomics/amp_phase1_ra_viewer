theme_clean <- function(base_size = 20, base_family = "") {
  theme_classic(
    base_size = base_size, base_family = base_family
  ) %+replace%
    theme(
      # axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.25),
      strip.background = element_rect(size = 0),
      panel.background = element_rect(size = 0.5, color = "black"),
      axis.line = element_blank()
    )
}
