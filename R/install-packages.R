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

# Load packages, and install them if they are not installed.
if (!require(pacman)) { install.packages("pacman") }
pacman::p_load(char = cran_packages)

github_packages <- c(
  "AnalytixWare/ShinySky",
  "baptiste/egg",
  "eclarke/ggbeeswarm",
  "mojaveazure/loomR"
)

for (github_package in github_packages) {
  devtools::install_github(github_package)
}

