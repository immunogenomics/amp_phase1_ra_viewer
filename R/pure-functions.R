scale_rows <- function(...) t(scale(t(...)))

extreme_n <- function(x, n = 10) {
  x <- sort(x)
  i <- floor(n / 2)
  j <- length(x) - i
  c(x[1:i], x[j:length(x)])
}
