# ampviewer
# Kamil Slowikowski
# 2018-01-26
#
# Fan Zhang Updated (added more visualization fuctions)
# 2018-01-29
# 
# This is an app for viewing single-cell RNA-seq data from AMP.
# Fan Zhang produced the input data (tSNE, clusters).
#
# getwd()
# setwd("/Users/fanzhang/Documents/HMS/amp/results/2018_01_26_ampviewer")

# Libraries -------------------------------------------------------------------

library(shiny)
library(shinysky) # devtools::install_github("AnalytixWare/ShinySky")

library(ggplot2)
library(ggbeeswarm)
library(magrittr)
library(dplyr)
# library(scales)
# library(RColorBrewer)
# library(pryr)
# library(egg)
# library(loomR) # devtools::install_github("mojaveazure/loomR")

# library(viridis)
# library(grid)
# library(Matrix)
# library(parallel)
# library(formattable)
# library(pheatmap)
# library(d3heatmap)

source("R/plot-tsne.R")
source("R/plot-box.R")
source("R/load-data.R")

# User interface --------------------------------------------------------------

# Call this function with all the regular navbarPage() parameters, plus a text parameter,
# if you want to add text to the navbar
navbarPageWithText <- function(..., text) {
  navbar <- navbarPage(...)
  textEl <- tags$div(class = "navbar-text pull-right", text)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]],
    textEl
  )
  navbar
}
  
ui <- fluidPage(

  tags$head(
    includeScript("google-analytics.js"),
    tags$link(
      rel = "stylesheet", type = "text/css", href = "app.css"
    ),
    tags$style("#tnse_marker_plot{min-width:500px;max-height:500px;}")
  ),

  navbarPageWithText(
    "AMP Phase I",
    source(file.path("R", "ui-tab-ra.R"), local = TRUE)$value,
    source(file.path("R", "ui-tab-about.R"), local = TRUE)$value,
    text = textOutput("mem_used")
  )

)

#

# Server ----------------------------------------------------------------------

# Debug
# input <- list(cell_type = "fibro", one_gene_symbol = "IFNB1")

server <- function(input, output) {
  
  output$tnse_marker_plot <- renderPlot({
    marker <- ifelse(
      input$one_gene_symbol != "" && input$one_gene_symbol %in% gene_symbols,
      input$one_gene_symbol,
      one_gene_symbol_default
    )
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- lf$matrix[,gene_ix]
    cell_ix <- which(cell_types == input$cell_type)
    plot_tsne(meta[cell_ix,], marker)
  })
  
  # output$marker_heatmap <- renderD3heatmap({
  #   markers <- get(sprintf("markers_%s", input$cell_type))
  #   d3heatmap(
  #     x               = -log10(markers),
  #     colors          = "Greys",
  #     yaxis_font_size = "14px"
  #   )
  # })
  
  # output$box_marker_plot_single <- renderPlot({
  #   marker <- ifelse(
  #     input$one_gene_symbol != "" && input$one_gene_symbol %in% gene_symbols,
  #     input$one_gene_symbol,
  #     one_gene_symbol_default
  #   )
  #   gene_ix <- which(gene_symbols == marker)
  #   meta$marker <- lf$matrix[,gene_ix]
  #   cell_ix <- which(cell_types == input$cell_type)
  #   plot_box(meta[cell_ix,], marker)
  # })
  
  output$box_marker_plot_all <- renderPlot({
    marker <- ifelse(
      input$one_gene_symbol != "" && input$one_gene_symbol %in% gene_symbols,
      input$one_gene_symbol,
      one_gene_symbol_default
    )
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- lf$matrix[,gene_ix]
    plot_box(meta, marker)
  })
  
  output$mem_used <- renderText({
       gdata::humanReadable(pryr::mem_used(), standard = "SI")
    })
  
  output$table <- renderDataTable(cluster_markers)
  
}

#

# Launch the app --------------------------------------------------------------

shinyApp(ui = ui, server = server)
