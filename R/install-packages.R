# This script installs the packages required for ampviewer.

installed_packages <- rownames(installed.packages())

cran_packages <- c(
  "RColorBrewer",
  "dplyr",
  "ggplot2",
  "magrittr",
  "pryr",
  "scales",
  "shiny"
)

# Load packages, and install them if they are not installed.
if (!"pacman" %in% installed_packages) { install.packages("pacman") }
pacman::p_load(char = cran_packages, install = TRUE, update = TRUE)

github_packages <- c(
  "AnalytixWare/shinysky",
  "baptiste/egg",
  "eclarke/ggbeeswarm",
  "mojaveazure/loomR"
)

for (github_package in github_packages) {
  p <- strsplit(github_package, "/")[[1]][2]
  if (!p %in% installed_packages) {
    devtools::install_github(github_package)
  }
}

