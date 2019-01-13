scale_rows <- function(...) t(scale(t(...)))

extreme_n <- function(x, n = 10) {
  x <- sort(x)
  i <- floor(n / 2)
  j <- length(x) - i
  c(x[1:i], x[j:length(x)])
}

slug <- function(x) {
  tolower(stringr::str_replace_all(x, " ", "_"))
}

#' Get the quantile breaks in a numeric vector.
#' @param x A numeric vector.
#' @param n The number of breaks.
#' @return A vector with unique breaks.
quantile_breaks <- function(x, n = 10) {
  breaks <- quantile(x, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

which_numeric_cols <- function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}

