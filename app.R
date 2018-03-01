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
# setwd("/Users/fanzhang/Documents/GitHub/ampviewer")

# https://github.com/timelyportfolio/d3treeR/issues/19#issuecomment-268110274
try({dev.off()})
pdf(NULL)

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
library(egg)
library(Matrix)
# library(loomR) # devtools::install_github("mojaveazure/loomR")

# library(viridis)
# library(grid)
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
    includeHTML("google-analytics.html"),
    tags$link(
      rel = "stylesheet", type = "text/css", href = "app.css"
    )
    #tags$style("#tnse_marker_plot{min-width:500px;max-height:500px;}")
  ),

  navbarPageWithText(
    "AMP Phase I",
    source(file.path("R", "ui-tab-ra.R"), local = TRUE)$value,
    source(file.path("R", "ui-tab-about.R"), local = TRUE)$value,
    text = uiOutput("navbar_right")
  )

)

#

# Server ----------------------------------------------------------------------

which_numeric_cols <- function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}

#' Test if a program is installed.
is_installed <- function(command) {
  !startsWith(system(
    command = sprintf("command -v %s || echo FALSE", command),
    intern = TRUE
  ), "FALSE")
}

#' Call pngquant to optimize a PNG file.
optimize_png <- function(filename) {
  if (is_installed("pngquant")) {
    opt_filename <- sprintf(
      "%s-fs8.png", substr(filename, 1, nchar(filename) - 4)
    )
    command <- sprintf(
      "pngquant --speed=1 --ext -fs8.png -- %s && mv -f %s %s",
      filename, opt_filename, filename
    )
    system(command)
  }
}

server <- function(input, output, session) {
  
  # Debug
  # input <- list(cell_type = "fibro", one_gene_symbol = "IFNB1")
  
  output$tnse_marker_plot <- renderImage({
    marker <- one_gene_symbol_default
    this_gene <- dg$gene[input$dg_table_rows_selected]
    if (length(this_gene) > 0) {
      marker <- this_gene
    }
    stopifnot(marker %in% gene_symbols)
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- lf$matrix[,gene_ix]
    
    stopifnot(input$cell_type %in% possible_cell_types)
    if (input$cell_type == "all") {
      cell_ix <- seq(nrow(meta))
      tsne_x <- "T1_all"
      tsne_y <- "T2_all"
    } else {
      cell_ix <- which(cell_types == input$cell_type)
      tsne_x <- "T1"
      tsne_y <- "T2"
    }
    
    temp_dir <- file.path("/tmp/ampviewer")
    dir.create(temp_dir, showWarnings = FALSE)
    outfile <- file.path(
      temp_dir,
      sprintf("tsne_%s_%s.png", input$cell_type, marker)
    )
    
    if (!file.exists(outfile) || file_test("-nt", "app.R", outfile)) {
      p <- plot_tsne(meta[cell_ix,], tsne_x, tsne_y, marker)
      ggsave(filename = outfile, plot = p, width = 10, height = 6, dpi = 100)
      # png(outfile, width = 600, height = 380)
      # print(p)
      # dev.off()
      optimize_png(outfile)
    }
    
    list(
      style = 'height: 100%; width: 100%; object-fit: contain',
      src = outfile,
      alt = sprintf("%s %s", input$cell_type, marker)
    )
  }, deleteFile = FALSE)
  
  # output$marker_heatmap <- renderD3heatmap({
  #   markers <- get(sprintf("markers_%s", input$cell_type))
  #   d3heatmap(
  #     x               = -log10(markers),
  #     colors          = "Greys",
  #     yaxis_font_size = "14px"
  #   )
  # })
  
  output$box_marker_plot_all <- renderImage({
    marker <- one_gene_symbol_default
    this_gene <- dg$gene[input$dg_table_rows_selected]
    if (length(this_gene) > 0) {
      marker <- this_gene
    }
    stopifnot(marker %in% gene_symbols)
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- lf$matrix[,gene_ix]
   
    temp_dir <- file.path("/tmp/ampviewer")
    dir.create(temp_dir, showWarnings = FALSE)
    outfile <- file.path(
      temp_dir,
      sprintf("beeswarm_%s.png", marker)
    )
    
    if (!file.exists(outfile) || file_test("-nt", "app.R", outfile)) {
      p <- plot_box(meta, marker)
      ggsave(filename = outfile, plot = p, width = 6, height = 9, dpi = 100)
      # png(outfile, width = 385, height = 520)
      # print(p)
      # dev.off()
      optimize_png(outfile)
    }
    
    list(
      style = 'height: 100%; width: 100%; object-fit: contain',
      src = outfile,
      alt = marker
    )
  }, deleteFile = FALSE)
  
  output$navbar_right <- renderUI({
    element <- ""
    hostname <- session$clientData$url_hostname
    if (any(startsWith(hostname, c("test.", "localhost", "127.0.0.1")))) {
      mem_used <- gdata::humanReadable(pryr::mem_used(), standard = "SI")
      git_hash <- system(
        command = "git log | head -n1 | cut -f2 -d' ' | cut -c1-7",
        intern = TRUE
      )
      element <- tags$div(
        tags$a(
          git_hash,
          href = sprintf(
            "https://github.com/immunogenomics/ampviewer/commit/%s", git_hash
          )
        ),
        mem_used
      )
    }
    element
  })
  
  # output$cluster_table <- renderTable(cluster_table)
  
  output$table <- renderDataTable(cluster_markers)
  
  output$dg_table <- DT::renderDataTable({
    numeric_cols <- colnames(dg)[which_numeric_cols(dg)]
    # Javascript-enabled table.
    DT::datatable(
      data = dg,
      colnames = c("Gene", "AUC", "-log10(P)", "Cluster"),
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = "Scroller",
      options = list(
        # columnDefs has bugs https://github.com/rstudio/DT/issues/311
        # columnDefs = list(list(targets = c(1), searchable = TRUE)),
        deferRender = TRUE,
        scrollY = 350, # Height in pixels
        scroller = TRUE,
        lengthMenu = FALSE,
        autoWidth = FALSE
      )
    ) %>%
      DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
}

# Launch the app --------------------------------------------------------------

shinyApp(ui = ui, server = server)


