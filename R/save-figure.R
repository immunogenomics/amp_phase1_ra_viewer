save_figure <- function(
  filename = NULL, width = 6, height = 5, dpi = 100,
  html_style = "height: 100%; width: 100%; object-fit: contain",
  html_alt = NULL,
  ggplot_function = NULL
) {
  out_dir <- "www/figures"
  dir.create(out_dir, showWarnings = FALSE)
  filename <- file.path(out_dir, filename)
  if (!file.exists(filename) || file_test("-nt", "app.R", filename)) {
    if (is.null(ggplot_function)) {
      stop("ggplot_function is NULL")
    }
    p <- ggplot_function()
    ggsave(
      filename = filename, plot = p,
      width = width, height = height, dpi = dpi
    )
    optimize_png(filename)
  }
  glue(
    '<img style="{style}" src="{src}" alt="{alt}"></img>',
    style = html_style,
    src = str_remove(filename, "www/"),
    alt = html_alt
  )
}