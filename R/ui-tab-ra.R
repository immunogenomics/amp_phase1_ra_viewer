tabPanel(
  "Rheumatoid Arthritis",
  
  tabsetPanel(
    
    tabPanel(
      h4("Search genes"),
    
      # Sidebar with a slider input for number of bins
      sidebarLayout(
        
        sidebarPanel(
          
          radioButtons(
            inputId  = "cell_type",
            label    = "Cell type:",
            choices  = possible_cell_types,
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
            p("Explore gene expresson in single-cell RNA-seq clusters."),
  
            # h4("The tSNE plot shows the expression of selected gene in the cells from the selected cell type."),
            plotOutput("tnse_marker_plot", height = "700px"),
            br(),
            # hr(),
            # h2("Wilcox -Log10 P for 1 vs all"),
            # p("Top 100 genes for each cluster. May be high or low in a",
            #  " cluster."),
            # d3heatmapOutput("marker_heatmap", height = "1600px")
            hr(),
            
            # h4("The expression of selected gene in the selected Cell Type subsets (in violinplot):"),
            # plotOutput("box_marker_plot_single", height = "500px"),
            # br(),
            # hr(),
            
            # h4("The expression of selected gene across all the subsets (in violinplot):"),
            plotOutput("box_marker_plot_all", height = "500px"),
            br()
            
          )
        )
        
      ) # sidebarLayout
      
    ), # tabPanel
    
    tabPanel(
      h4("Marker genes"),
      mainPanel(
      br(),
      h4("Table of identified subsets marker genes (19 subsets in all)"),
      br(),
      fluidRow(
        column(12,
               dataTableOutput('table')
        )
      )
      )
    )
    # tabPanel(
    #   "Boxplot",
    #   mainPanel(
    #     fluidRow(
    #       plotOutput("box_marker_plot", height = "800px")
    #     )
    #   )
    # ) # tabPanel
    
    
  ) # tabsetPanel
  
) # tabPanel
