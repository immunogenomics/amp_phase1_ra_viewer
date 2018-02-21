tabPanel(
  "Rheumatoid Arthritis",
  
  tabsetPanel(
    
    tabPanel(
      h4("Search genes"),
        
        # Show a plot of the generated distribution
        fluidPage(
          fluidRow(
  
            column(width = 9, plotOutput("tnse_marker_plot", height = "400px")),
            
            column(width = 3,
              wellPanel(radioButtons(
                inputId  = "cell_type",
                # inline   = TRUE,
                label    = "Cell type:",
                choices  = possible_cell_types,
                selected = "fibro"
              ))
            )
            
          ),
          hr(),
          fluidRow(
            
            column(width = 6, plotOutput("box_marker_plot_all", height = "500px")),
            column(width = 6, DT::dataTableOutput("dg_table"))
            
          )
        )
      
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
