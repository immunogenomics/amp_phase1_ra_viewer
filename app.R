# ampviewer
# Kamil Slowikowski
# 2018-01-26
#
# Fan Zhang Updated (added more visualization fuctions)
# 2018-01-29
# 
# This is a minimal app for viewing single-cell RNA-seq data from AMP.
# Fan Zhang produced the input data (tSNE, clusters).

# devtools::install_github("AnalytixWare/ShinySky")
# getwd()
# setwd("/Users/fanzhang/Documents/HMS/amp/results/2018_01_26_ampviewer")

# Libraries -------------------------------------------------------------------

library(shiny)
library(shinysky)

library(grid)
library(ggplot2)
library(scales)
library(viridis)
library(egg)

library(Matrix)
# library(parallel)
# library(formattable)

#library(pheatmap)
# library(d3heatmap)

meta_colors <- list(
  "fine_cluster" = c(
    "CF1" = "#6BAED6",
    "CF2" = "#08306B",
    "CF3" = "#DEEBF7",
    "CF4" = "grey",
    "CT1" = "#FEB24C",
    "CT2" = "#8C510A",
    "CT3" = "brown",
    "CT4" = "#FFFF33",
    "CT5" = "#C7EAE5",
    "CT6" = "#003C30",
    "CT7" = "#35978F",
    "CB1" = "#FCBBA1",
    "CB2" = "#CB181D", #FB6A4A #A50F15
    "CB3" = "#67000D",
    "CB4" = "#FB9A99",
    "CM1" = "#AE017E",
    "CM2" = "#F768A1",
    "CM3" = "#FDE0EF", #FCC5C0
    "CM4" = "#49006A"
  )
)

# Prepare data ----------------------------------------------------------------

get_markers <- function(log2cpm, clusters) {
  # Get a Wilcox p-value for each gene and each cluster. 
  retval <- parallel::mclapply(X = unique(clusters), FUN = function(cluster) {
     ix_x <- which(clusters == cluster)
     ix_y <- which(clusters != cluster)
     pvals <- apply(log2cpm, 1, function(row) {
      wilcox.test(
        x = row[ix_x],
        y = row[ix_y],
        alternative = "two.sided"
      )$p.value
    })
  }, mc.cores = 4)
  retval <- do.call(cbind, retval)
  # Take the top 100 genes for each cluster.
  x <- apply(retval, 2, function(x) x < sort(x)[100])
  # Choose the genes that are specific to 1 cluster.
  best_markers <- retval[which(rowSums(x) == 1),]
  colnames(best_markers) <- unique(clusters)
  best_markers
}

cell_types <- c(
  "B cell"     = "bcell",
  "T cell"     = "tcell",
  "Monocyte"   = "mono",
  "Fibroblast" = "fibro"
)

# Read 4 datasets: bcell, tcell, mono, fibro
# Preprocess into a file for quick loading.
data_file <- "data/shiny.rda"
if (file.exists(data_file)) {
  load(data_file)
} else {
  e <- environment()
  for (cell_type in cell_types) {
    log2cpm <- readRDS(sprintf("data/%s_exp.rds", cell_type))
    log2cpm <- Matrix(log2cpm)
    meta    <- readRDS(sprintf("data/%s_sc_label.rds", cell_type))
    nonzero <- rownames(log2cpm)[
       which(rowSums(log2cpm > 0) > 10)
    ]
    # log2cpm_filter <- log2cpm[nonzero,]
    # markers <- get_markers(log2cpm_filter, meta$cluster)
    assign(sprintf("log2cpm_%s", cell_type), log2cpm, envir = e)
    assign(sprintf("meta_%s", cell_type), meta, envir = e)
    assign(sprintf("nonzero_%s", cell_type), nonzero, envir = e)
    # assign(sprintf("markers_%s", cell_type), markers, envir = e)
  }
  #all_cell_types <- cbind.data.frame(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
  all_meta <- rbind.data.frame(meta_fibro, meta_bcell, meta_tcell, meta_mono)
  rm(log2cpm)
  rm(meta)
  rm(nonzero)
  # rm(markers)
  gene_symbols <- unique(c(
    nonzero_bcell, nonzero_fibro, nonzero_mono, nonzero_tcell
  ))
  save.image(data_file)
}

# all_cell_types <- cbind.data.frame(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
# all_meta <- rbind.data.frame(meta_fibro, meta_bcell, meta_tcell, meta_mono)
# all(colnames(all_cell_types) == rownames(all_meta))

one_gene_symbol_default <- "HLA-DRA"

#

# Functions -------------------------------------------------------------------

quantile_breaks <- function(x, n = 10) {
  breaks <- quantile(x, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

plot_tsne <- function(log2cpm, dat, marker) {
  dat$marker <- as.numeric(log2cpm[marker,])
  n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 3.5
  fill_values <- quantile_breaks(dat$marker, n = 9)
  fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  theme_tsne <- theme_bw(base_size = 22) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(size = 25,  face="bold")
  )
  # Put the name of the marker gene in the upper left corner
  dat_text <- data.frame(
    x     = -Inf,
    y     = Inf,
    label = marker
  )
  p1 <- ggplot() +
    geom_point(
      data    = dat[order(dat$marker),],
      mapping = aes(x = T1, y = T2, fill = marker),
      size    = point_size,
      shape   = 21,
      stroke  = 0.15
    ) +
    # geom_text(
    #   data    = dat_text,
    #   mapping = aes(x, y, label = label),
    #   size    = 9,
    #   hjust   = -0.05,
    #   vjust   = 1.25
    # ) +
    scale_fill_gradientn(
      # Linear scale
      # colours = fill_palette,
      # Quantile scale
      colours = colorRampPalette(fill_palette)(length(fill_values)),
      values  = fill_values,
      breaks  = scales::pretty_breaks(n = 4),
      name    = bquote("Log"[2]~"(CPM+1)  ")
    ) +
    guides(
      fill  = guide_colorbar(barwidth = 10, barheight = 1),
      alpha = "none"
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(marker) +
    theme_tsne
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  p2 <- ggplot() +
    geom_point(
      data    = dat[sample(nrow(dat)),],
      mapping = aes(x = T1, y = T2, fill = cluster),
      size    = point_size,
      shape   = 21,
      stroke  = 0.15
    ) +
    # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
    labs(x = NULL, y = NULL) +
    ggtitle("Identified clusters") +
    theme_tsne
  bottom_text <- sprintf(
    "%s cells, %s (%s%%) nonzero cells",
    nrow(dat),
    n_nonzero,
    signif(100 * n_nonzero / nrow(dat), 3)
  )
  egg::ggarrange(
    bottom = textGrob(
      label = bottom_text, gp = gpar(fontsize = 20)
    ),
    plots = list(p1, p2), ncol = 2
  )
}
# plot_tsne(log2cpm_fibro, meta_fibro, "CD3D")

plot_box <- function(log2cpm_marker, dat, marker) {
  #dat$marker <- as.numeric(log2cpm[marker,])
  dat$marker <- as.numeric(log2cpm_marker)
  theme_box <- theme_bw(base_size = 22) + theme(
    legend.position = "bottom",
    # axis.text       = element_blank(),
    # axis.ticks      = element_blank(),
    # panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(size = 25,  face="bold")
  )
  dat$cluster <- factor(dat$cluster)
  ggplot(
    data=dat, 
    aes(x=cluster, 
        y=marker, 
        fill=cluster)) +
    # geom_boxplot() +
    geom_violin() +
    geom_jitter(height = 0, width = 0.25, color = "dimgrey", size = 0.5) + 
    labs(
      x = NULL,
      y    = bquote("Log"[2]~"(CPM+1)  "),
      title = marker
      # subtitle = tsne_subtitle
    ) +
    # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
    theme_box
  
  # proportion <- rep(0, length(table(dat$cluster)))
  # for (i in 1:length(table(dat$cluster))){
  #   proportion[i] <- sum(dat$cluster == i & dat$marker > 0)/ (table(dat$cluster)[i])
  #   # print(sum(dat$cluster == i & dat$marker > 0))
  #   # print(table(dat$cluster)[i])
  # }
  # dat_pro <- data.frame(
  #   cluster = as.character(seq(1, length(table(dat$cluster)))),
  #   nonzero = proportion
  # )
  # p2 <- ggplot(
  #   data=dat_pro, 
  #   aes(x=cluster, y= nonzero, fill = cluster)
  #   # aes(x=cluster, y= percent(nonzero), fill = cluster)
  #   ) +
  #   geom_bar(stat="identity", position = "stack") +
  #   labs(
  #     x = NULL,
  #     y = "Percent nonzero",
  #     title = marker
  #   ) +
  #   scale_y_continuous(labels = percent) +
  #   # scale_fill_brewer(type = "qual", palette = "Set3", name = "Cluster") +
  #   scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  #   theme_box
  # p3 <- ggplot()
  # bottom_text <- sprintf(
  #   "%s is expressing the highest percent of nonzero cells in cluster %s.",
  #   marker,
  #   as.integer(dat_pro$cluster[which(dat_pro$nonzero == max(dat_pro$nonzero))])
  # )
  # egg::ggarrange(
  #   # bottom = textGrob(
  #   #   label = bottom_text, 
  #   #   gp = gpar(fontsize = 20, fontface="bold")
  #   # ),
  #   plots = list(p1, p2), ncol = 2, widths = 2.7:1, height = 2.2:1,
  #   hr()
  # )
}

# For testing, it's nice to have a little snippet here.
# plot_tsne(
#   log2cpm = log2cpm_mono,
#   dat     = meta_mono,
#   marker  = one_gene_symbol_default
# )

#

# User interface --------------------------------------------------------------

# Call this function with all the regular navbarPage() parameters, plus a text parameter,
# if you want to add text to the navbar
navbarPageWithText <- function(..., text) {
  navbar <- navbarPage(...)
  textEl <- tags$div(class = "navbar-text pull-right", text)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]], textEl)
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
  
  # Application title
  navbarPageWithText(
    "AMP Phase I",
    
    tabPanel(
      "Rheumatoid Arthritis",
      
      # h3("Integration of Single-cell Transcriptomic and Proteomic 
      #    Immune Profiling Identifies Pathogenic Pathways in Rheumatoid Arthritis"),
      
      p("Explore gene expression in single-cell RNA-seq clusters."),
      br(),
      
      tabsetPanel(
        
        tabPanel(
          "Search genes",
        
          # Sidebar with a slider input for number of bins
          sidebarLayout(
            
            sidebarPanel(
              
              # selectInput(
              #   inputId  = "cell_type",
              #   label    = "Cell type:",
              #   choices  = cell_types,
              #   selected = "fibro"
              # ),
              radioButtons(
                inputId = "cell_type",
                label   = "Cell type:",
                choices = cell_types,
                selected = "fibro"
              ),
              
              strong("Gene:"),
              
              textInput.typeahead(
                id = "one_gene_symbol",
                placeholder = one_gene_symbol_default,
                local = data.frame(name = gene_symbols),
                limit = 10,
                valueKey = "name",
                tokens = seq_along(gene_symbols),
                template = HTML("<p>{{name}}</p>")
              ),
              
              # Sidebar is fluid and 3/12 units wide.
              width = 3,
              
              hr(),
              h5("Contact"),
              p(
                "This site is made by", 
                a("Kamil Slowikowski", href = "mailto:kslowikowski@fas.harvard.edu"),
                "and",
                a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
                "Please contact",
                a("Fan", href = "mailto:fanzhang@broadinstitute.org"),
                "if you have any questions, requests, or comments on the analysis and results."
              )
            ),
            
            # Show a plot of the generated distribution
            mainPanel(
              fluidRow(
                h4("The tSNE plot shows the expression of selected gene in the cells from the selected cell type."),
                plotOutput("tnse_marker_plot", height = "700px"),
                br(),
                # hr(),
                # h2("Wilcox -Log10 P for 1 vs all"),
                # p("Top 100 genes for each cluster. May be high or low in a",
                #  " cluster."),
                # d3heatmapOutput("marker_heatmap", height = "1600px")
                hr(),
                h4("The expression of selected gene in the selected Cell Type subsets (in violinplot):"),
                plotOutput("box_marker_plot_single", height = "500px"),
                br(),
                hr(),
                h4("The expression of selected gene across all the subsets (in violinplot):"),
                plotOutput("box_marker_plot_all", height = "500px"),
                br()
                
              )
            )
            
          ) # sidebarLayout
          
        ) # tabPanel
        
        # tabPanel(
        #   "Boxplot",
        #   mainPanel(
        #     fluidRow(
        #       plotOutput("box_marker_plot", height = "800px")
        #     )
        #   )
        # ) # tabPanel
        
        
      ) # tabsetPanel
      
    ),
    
    tabPanel(
      
      "About",
      
      mainPanel(
        h1("Accelerating Medicines Partnerships (AMP)"),
        p(
          "The",
          a("Accelerating Medicines Partnership (AMP)",
            href = "https://www.nih.gov/research-training/accelerating-medicines-partnership-amp"),
          " is a public-private partnership between the National Institutes of",
          " Health (NIH), the U.S. Food and Drug Administration (FDA), 10",
          " biopharmaceutical companies and multiple non-profit organizations",
          " to transform the current model for developing new diagnostics and",
          " treatments by jointly identifying and validating promising",
          " biological targets for therapeutics. The ultimate goal is to",
          " increase the number of new diagnostics and therapies for patients",
          " and reduce the time and cost of developing them."
        ),
        h2("Disclaimer"),
        p(
          "Data presented on this page is from Phase 1 of the AMP partership.",
          " Currently, this is private data meant to be shared internally,",
          " only with consortium members."
        ),
        p(
          strong(
            "Sharing any data from this site with anyone outside of the",
            " AMP partnership is prohibited."
          )
        ),
        p(
          "This website is an experiment in providing early access to",
          " preliminary data analysis results. The content of this site is",
          " subject to change at any time without notice. We hope that you",
          " find it useful, but we provide it 'as is' without warranty of",
          " any kind, express or implied."
        ),
        h2("Contact"),
        p(
          "This site is made by", 
                    a("Kamil Slowikowski", href = "mailto:kslowikowski@fas.harvard.edu"),
          "and",
                    a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
          "Please contact",
                    a("Fan", href = "mailto:fanzhang@broadinstitute.org"),
          "if you have any questions, requests, or comments on the analysis and results."
        )
        
      ) # mainPanel
      
    ), # tabPanel
    
    text = textOutput("mem_used")
  ) # navbarPage
  
) # fluidPage

#

# Server ----------------------------------------------------------------------

server <- function(input, output) {
  
  output$tnse_marker_plot <- renderPlot({
    log2cpm <- get(sprintf("log2cpm_%s", input$cell_type))
    dat     <- get(sprintf("meta_%s", input$cell_type))
    marker  <- ifelse(
      input$one_gene_symbol != "",
      input$one_gene_symbol,
      one_gene_symbol_default
    )
    # Don't allow selecting genes that are not present in the data.
    if (! marker %in% rownames(log2cpm)) {
      marker <- rownames(log2cpm)[1]
    }
    plot_tsne(
      log2cpm = log2cpm,
      dat     = dat,
      marker  = marker
    )
  })
  
  # output$marker_heatmap <- renderD3heatmap({
  #   markers <- get(sprintf("markers_%s", input$cell_type))
  #   d3heatmap(
  #     x               = -log10(markers),
  #     colors          = "Greys",
  #     yaxis_font_size = "14px"
  #   )
  # })
  
  output$box_marker_plot_single <- renderPlot({
    log2cpm <- get(sprintf("log2cpm_%s", input$cell_type))
    dat     <- get(sprintf("meta_%s", input$cell_type))
    marker  <- ifelse(
      input$one_gene_symbol != "",
      input$one_gene_symbol,
      one_gene_symbol_default
    )
    # Don't allow selecting genes that are not present in the data.
    if (! marker %in% rownames(log2cpm)) {
      marker <- rownames(log2cpm)[1]
    }
    plot_box(
      log2cpm_marker = as.numeric(log2cpm[marker,]),
      dat     = dat,
      marker  = marker
    )
  })
  
  output$box_marker_plot_all <- renderPlot({
    marker  <- ifelse(
      input$one_gene_symbol != "",
      input$one_gene_symbol,
      one_gene_symbol_default
    )
    # Don't allow selecting genes that are not present in the data.
    if (! marker %in% rownames(log2cpm_fibro)) {
      marker <- rownames(log2cpm_fibro)[1]
    }
    log2cpm_marker <- c(
      as.numeric(log2cpm_fibro[marker,]),
      as.numeric(log2cpm_bcell[marker,]),
      as.numeric(log2cpm_tcell[marker,]),
      as.numeric(log2cpm_mono[marker,])
    )
    #log2cpm <- get(sprintf("all_cell_types"))
    dat <- get(sprintf("all_meta"))
    plot_box(
      log2cpm_marker = log2cpm_marker,
      #log2cpm = log2cpm,
      dat     = dat,
      marker  = marker
    )
  })
  
  output$mem_used <- renderText({
    gdata::humanReadable(pryr::mem_used(), standard = "SI")
  })
  
}

#

# Launch the app --------------------------------------------------------------

shinyApp(ui = ui, server = server)

