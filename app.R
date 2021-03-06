# Data viewer for AMP Phase I RA
# creator: Kamil Slowikowski
# contributors: Fan Zhang, Joseph Mears
# 
# https://github.com/timelyportfolio/d3treeR/issues/19#issuecomment-268110274
try({dev.off()})
pdf(NULL)

# Should we save each figure and reuse it instead of regenerating each time?
USE_CACHE <- TRUE

# Dependencies ----------------------------------------------------------------

source("R/install-packages.R")
source("R/meta-colors.R")
source("R/optimize-png.R") 
source("R/save-figure.R")
source("R/theme-clean.R")
source("R/pure-functions.R")
source("R/plot-tsne.R")
source("R/plot-tsne-cytof.R")
source("R/plot-box.R")
source("R/plot-bulk-dots.R")
source("R/plot-bulk-single-cca.R")
source("R/plot-cca-scores.R")
source("R/load-data.R")

# User interface --------------------------------------------------------------

# Call this function with all the regular navbarPage() parameters, plus a text parameter,
# if you want to add text to the navbar
navbarPageWithText <- function(title, ..., id = NULL, text) {
  navbar <- navbarPage(title = title, ..., id = id)
  textEl <- tags$div(class = "navbar-text pull-right", text)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]],
    textEl
  )
  navbar
}
  
ui <- function(request) {
fluidPage(

  tags$head(
    includeHTML("google-analytics.html"),
    tags$link(
      rel = "stylesheet", type = "text/css", href = "app.css"
    ),
    tags$script(src = "utils.js")
    #tags$style("#tnse_marker_plot{min-width:500px;max-height:500px;}")
  ),

  navbarPageWithText(
    "AMP RA Phase I",
    source(file.path("R", "ui-tab-about.R"), local = TRUE)$value,
    source(file.path("R", "ui-tab-ra.R"), local = TRUE)$value,
    text = uiOutput("navbar_right")
  ),
  
  #bookmarkButton(),
  
  HTML(
    "<footer class='myfooter page-footer'>
    <div class='text-center'>
    This website was created by <a href='https://slowkow.com'>Kamil Slowikowski</a>,
    with help from <a href='https://immunogenomics.hms.harvard.edu/people/fan-zhang'>Fan Zhang</a>
    and <a href='https://immunogenomics.hms.harvard.edu/people/joseph-mears'>Joseph Mears</a>
    </div>
    </footer>
    "
  )

)
}

#

# Server ----------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Debug
  # input <- list(
  #   cell_type = "Fibroblast",
  #   one_gene_symbol = "IFNB1",
  #   dg_table_rows_selected = "IFNB1",
  #   one_protein_symbol = "CD90",
  #   cytof_table_rows_selected = "CD90",
  #   bulk_single_cca_xaxis = 1,
  #   bulk_single_cca_yaxis = 2
  # )

  output$tnse_marker_plot <- renderText({
    marker <- one_gene_symbol_default
    this_gene <- as.character(dg$gene[input$dg_table_rows_selected])
    if (length(this_gene) > 0) {
      marker <- this_gene
    }
    stopifnot(marker %in% gene_symbols)
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- as.numeric(log2cpm[gene_ix,])
    stopifnot(input$rnaseq_cell_type %in% possible_cell_types_rna)
    if (input$rnaseq_cell_type == "all") {
      cell_ix <- seq(nrow(meta))
      tsne_x <- "T1_all"
      tsne_y <- "T2_all"
    } else {
      cell_ix <- which(meta$cell_type == input$rnaseq_cell_type)
      tsne_x <- "T1" 
      tsne_y <- "T2"
    }
    save_figure(
      cache = USE_CACHE,
      filename = glue(
        "ampra1_scrnaseq_tsne_{celltype}_{marker}.jpg",
        celltype = slug(input$rnaseq_cell_type), marker = marker
      ),
      width = 14, height = 8, dpi = 100,
      html_alt = sprintf("%s %s", input$rnaseq_cell_type, marker),
      ggplot_function = function() {
        plot_tsne(meta[cell_ix,], tsne_x, tsne_y, title = marker)
      }
    )
  })
  
  # output$marker_heatmap <- renderD3heatmap({
  #   markers <- get(sprintf("markers_%s", input$rnaseq_cell_type))
  #   d3heatmap(
  #     x               = -log10(markers),
  #     colors          = "Greys",
  #     yaxis_font_size = "14px"
  #   )
  # })
  
  output$box_marker_plot_all <- renderText({
    marker <- one_gene_symbol_default 
    this_gene <- as.character(dg$gene[input$dg_table_rows_selected])
    if (length(this_gene) > 0) {
      marker <- this_gene
    }
    stopifnot(marker %in% gene_symbols)
    gene_ix <- which(gene_symbols == marker)
    meta$marker <- as.numeric(log2cpm[gene_ix,])
    save_figure(
      cache = USE_CACHE,
      filename = glue("ampra1_scrnaseq_bar_{marker}.jpg", marker = marker),
      width = 6, height = 9, dpi = 100,
      html_alt = marker,
      ggplot_function = function() { plot_box(meta, marker) }
    )
  })
  
  output$bulk_dots <- renderText({
    marker <- one_gene_symbol_default
    this_gene <- as.character(dg$gene[input$dg_table_rows_selected])
    if (length(this_gene) > 0) {
      marker <- this_gene
    }
    # stopifnot(marker %in% gene_symbols)
    retval <- glue(
      "<p class='alert alert-warning'><i>{gene}</i> is not present in the bulk RNA-seq data.</p>",
      gene = marker
    )
    if (marker %in% rownames(b_log2tpm)) {
      b_meta$marker <- as.numeric(b_log2tpm[marker,])
      retval <- save_figure(
        cache = USE_CACHE,
        filename = glue("ampra1_rnaseq_dots_{marker}.jpg", marker = marker),
        width = 9, height = 8, dpi = 100,
        html_alt = marker,
        ggplot_function = function() { plot_bulk_dots(b_meta, marker) }
      )
    }
    retval
  })
  
  # output$bulk_single_cca <- renderText({
  #   marker <- one_gene_symbol_default
  #   this_gene <- as.character(dg$gene[input$dg_table_rows_selected])
  #   if (length(this_gene) > 0) {
  #     marker <- this_gene
  #   }
  #   stopifnot(marker %in% gene_symbols)
  #   gene_ix <- which(gene_symbols == marker)
  #   meta$marker <- as.numeric(log2cpm[gene_ix,])
  #   sc_marker <- meta$marker
  #   names(sc_marker) <- as.character(meta$cell_name)
  #   
  #   dat_cca$marker <- c(
  #     as.numeric(b_log2tpm[marker,]),
  #     as.numeric(sc_marker[cca_bs_ynames])
  #   )
  #   dimx <- input$bulk_single_cca_xaxis
  #   dimy <- input$bulk_single_cca_yaxis
  #   if (!dimx %in% 1:10) {
  #     dimx <- 1
  #   }
  #   if (!dimy %in% 1:10) {
  #     dimy <- 2
  #   }
  #   save_figure(
  #     cache = USE_CACHE,
  #     filename = glue(
  #       "ampra1_cca_rnaseq_scrnaseq_cv{x}_cv{y}_{marker}.png",
  #       x = dimx,
  #       y = dimy,
  #       marker = marker
  #     ),
  #     width = 8, height = 8, dpi = 100,
  #     html_alt = marker,
  #     ggplot_function = function() {
  #       plot_bulk_single_cca(
  #         dat_cca = dat_cca,
  #         x = as.character(dimx),
  #         y = as.character(dimy),
  #         marker = marker
  #       )
  #     }
  #   )
  # })
  
  # output$bulk_single_cca_scores <- renderText({
  #   dimx <- input$bulk_single_cca_xaxis
  #   dimy <- input$bulk_single_cca_yaxis
  #   if (!dimx %in% 1:10) {
  #     dimx <- 1
  #   }
  #   if (!dimy %in% 1:10) {
  #     dimy <- 2
  #   }
  #   save_figure(
  #     cache = USE_CACHE,
  #     filename = glue(
  #       "ampra1_cca_rnaseq_scrnaseq_scores_{x}_{y}.png",
  #       x = dimx,
  #       y = dimy
  #     ),
  #     width = 8, height = 7, dpi = 100,
  #     html_alt = "CCA scores",
  #     ggplot_function = function() {
  #       plot_cca_scores(
  #         cca_bs, cv = as.integer(dimx), n = 20
  #       ) +
  #       plot_cca_scores(
  #         cca_bs, cv = as.integer(dimy), n = 20
  #       )
  #     }
  #   )
  # })
  
  output$navbar_right <- renderUI({
    element <- ""
    hostname <- session$clientData$url_hostname
    if (any(startsWith(hostname, c("test.", "localhost", "127.0.0.1")))) {
      mem_used <- gdata::humanReadable(pryr::mem_used(), standard = "SI")
      git_date <- system(
        command = "git log -1 --pretty=format:'%ci' | cut -f1 -d' '",
        intern = TRUE
      )
      git_hash <- system(
        command = "git log -1 --pretty=format:'%h'",
        intern = TRUE
      )
      element <- tags$div(
        git_date,
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
  
  
  output$table <- renderDataTable(cluster_markers)
  
  output$dg_table <- DT::renderDataTable({
    numeric_cols <- colnames(dg)[which_numeric_cols(dg)]
    # Javascript-enabled table.
    DT::datatable(
      data = dg,
      colnames = c("Gene", "Single-cell RNA-seq cluster", "AUC", "% nonzero"),
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
  
  
  
  output$tnse_cytof <- renderText({
    marker <- one_protein_symbol_default
    # this_protein <- as.character(cytof_summarize$protein[input$cytof_table_rows_selected])
    this_protein <- as.character(input$cytof_marker)
    if (length(this_protein) > 0) {
      marker <- this_protein
    }
    stopifnot(marker %in% protein_symbols)
    # stopifnot(input$cytof_cell_type %in% possible_cell_types_cytof)
    if (!input$cytof_cell_type %in% possible_cell_types_cytof) {
      # textOutput("The mass cytometry data were analyzed per cell type.")
      return("<p class='alert alert-danger'><strong>Error:</strong> Please select one cell type, not 'All cells'.</p>")
    }
    cytof_all$marker <- as.numeric(cytof_all[, which(colnames(cytof_all) == marker)])
    cell_ix <- which(cytof_all$cell_type == input$cytof_cell_type)
    tsne_x <- "SNE1" 
    tsne_y <- "SNE2"
      
    save_figure(
      cache = USE_CACHE,
      filename = glue(
        "ampra1_cytof_tsne_{celltype}_{marker}.jpg",
        celltype = slug(input$cytof_cell_type), marker = marker
      ),
      width = 14, height = 9, dpi = 100,
      html_alt = sprintf("%s %s", input$cytof_cell_type, marker),
      ggplot_function = function() {
        plot_tsne_cytof(cytof_all[cell_ix,], tsne_x, tsne_y, title = marker)
      }
    )
  })
  
  output$cytof_table <- DT::renderDataTable({
    numeric_cols <- names(which(sapply(cytof_summarize, is.numeric)))
    # Javascript-enabled table.
    DT::datatable(
      data       = cytof_summarize,
      colnames   = c("Cell Type", "Cluster", "Protein", "Percent Non-zero"),
      selection  = "single",
      rownames   = FALSE,
      filter     = list(position = "top", plain = TRUE),
      style      = "default",
      extensions = "Scroller",
      options    = list(
        # columnDefs has bugs https://github.com/rstudio/DT/issues/311
        # columnDefs = list(list(targets = c(1), searchable = TRUE)),
        deferRender = TRUE,
        scrollY     = 350, # Height in pixels
        scroller    = TRUE,
        lengthMenu  = FALSE,
        autoWidth   = FALSE
      )
    ) %>%
      DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
  
}

# Launch the app --------------------------------------------------------------

# Bookmarks would be great, but:
#
# The data table records all of the rows as "input" for some reason, so the
# URL encoded values have 100s or 1000s of rows from the table.
#
# enableBookmarking(store = "url")

shinyApp(ui = ui, server = server)


