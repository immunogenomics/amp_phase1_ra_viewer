# This script installs the packages required for ampviewer.

cran_packages <- c(
  "RColorBrewer",
  "dplyr",
  "ggplot2",
  "magrittr",
  "pryr",
  "scales",
  "shiny"
)

install.packages(cran_packages)

github_packages <- c(
  "AnalytixWare/ShinySky",
  "baptiste/egg",
  "eclarke/ggbeeswarm",
  "mojaveazure/loomR"
)

for (github_package in github_packages) {
  devtools::install_github(github_package)
}

