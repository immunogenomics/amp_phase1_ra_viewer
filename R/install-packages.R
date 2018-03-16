# This script installs the packages required for ampviewer.

installed_packages <- rownames(installed.packages())

cran_packages <- c(
  "RColorBrewer",
  "dplyr",
  "magrittr",
  "pryr",
  "hdf5r",
  "scales",
  "shiny",
  "DT",
  "lme4",
  "gdata"
)

# Load packages, and install them if they are not installed.
if (!"pacman" %in% installed_packages) { install.packages("pacman") }
pacman::p_load(char = cran_packages, install = TRUE, update = TRUE)

github_packages <- c(
  "AnalytixWare/shinysky",
  "eclarke/ggbeeswarm",
  "mojaveazure/loomR",
  "tidyverse/ggplot2",
  "thomasp85/patchwork"
)

for (github_package in github_packages) {
  p <- strsplit(github_package, "/")[[1]][2]
  if (!p %in% installed_packages) {
    devtools::install_github(github_package)
  }
}

