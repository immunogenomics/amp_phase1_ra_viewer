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
                selected = "all"
              ))
              # tableOutput("cluster_table")
            )
            
          ),
          hr(),
          fluidRow(
            column(
              width = 5,
              plotOutput("box_marker_plot_all", height = "500px")
            ),
            column(
              width = 7,
              DT::dataTableOutput("dg_table", height = "350px")
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 5,
              plotOutput("bulk_dots", height = "200px")
            ),
            column(
              width = 7,
              plotOutput("bulk_single_cca", height = "300px")
            )
          )
        )
      
    ), # tabPanel
    
    tabPanel(
      h4("Marker genes"),
      fluidPage(
        h4("Table of identified subsets marker genes (19 subsets in all)"),
        br(),
        fluidRow(
          column(12,
            dataTableOutput('table')
          )
        )
      )
    ),
    
    br()
    
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
