tabPanel(
  # "Rheumatoid Arthritis",
  "Data Viewer",

    fluidRow(
      column(
        width = 12,
        wellPanel(
          radioButtons(
          inputId  = "cell_type",
          inline   = TRUE,
          label    = "Cell type:",
          choices  = possible_cell_types_rna,
          selected = "Fibroblast"
       )
     )
    )
    ),
  
  
  tabsetPanel(selected = "Single-cell RNA-seq",
    
    tabPanel(
      "Single-cell RNA-seq",
        
        # Show a plot of the generated distribution
        fluidPage(
          fluidRow(
  
            column(
              width = 9,
              div(
                htmlOutput("tnse_marker_plot", height = "400px"),
                style = "height: 400px;"
              )
            )
            
            # column(
            #   width = 3,
            #   wellPanel(radioButtons(
            #     inputId  = "cell_type",
            #     # inline   = TRUE,
            #     label    = "Cell type:",
            #     choices  = possible_cell_types_rna,
            #     selected = "all"
            #   ))
              # tableOutput("cluster_table")
            )
            
          ),
          hr(),
          fluidRow(
            column(
              width = 5,
              htmlOutput("box_marker_plot_all", height = "500px")
            ),
            column(
              width = 7,
              DT::dataTableOutput("dg_table", height = "350px")
            )
          )
          # hr(),
          # fluidRow(
          #   column(
          #     width = 5,
          #     htmlOutput("bulk_dots", height = "500px")
          #   )
          # )
        #)
      
    ), # tabPanel
    
    # tabPanel(
    #   "CCA integration",
    #   fluidPage(
    #     fluidRow(
    #       # column(
    #       #   width = 5,
    #       #   fluidRow(
    #       #     column(width = 3, selectInput(
    #       #       inputId = "bulk_single_cca_xaxis",
    #       #       label = "x-axis",
    #       #       choices = 1:10,
    #       #       selected = 1
    #       #     )),
    #       #     column(width = 3, selectInput(
    #       #       inputId = "bulk_single_cca_yaxis",
    #       #       label = "y-axis",
    #       #       choices = 1:10,
    #       #       selected = 2
    #       #     ))
    #       #   )
    #       # ),
    #       
    #       column(
    #         width = 5,
    #         shinydashboard::box(
    #           width = 12,
    #           title = "Controls", 
    #           splitLayout(
    #             # selectInput(
    #             #   inputId = "bulk_single_cca_xaxis",
    #             #   label = "x-axis",
    #             #   choices = 1:10,
    #             #   selected = 1
    #             # ),
    #             # selectInput(
    #             #   inputId = "bulk_single_cca_yaxis",
    #             #   label = "y-axis",
    #             #   choices = 1:10,
    #             #   selected = 2
    #             # )
    #             # numericInput(
    #             #   inputId = "bulk_single_cca_xaxis",
    #             #   label = "x-axis",
    #             #   value = 1,
    #             #   min = 1,
    #             #   max = 10,
    #             #   step = 1
    #             # ),
    #             # numericInput(
    #             #   inputId = "bulk_single_cca_yaxis",
    #             #   label = "y-axis",
    #             #   value = 2,
    #             #   min = 1,
    #             #   max = 10,
    #             #   step = 1
    #             # )
    #             sliderInput(
    #               inputId = "bulk_single_cca_xaxis",
    #               label = "x-axis",
    #               value = 1,
    #               min = 1,
    #               max = 10,
    #               step = 1,
    #               round = TRUE,
    #               ticks = TRUE
    #             ),
    #             sliderInput(
    #               inputId = "bulk_single_cca_yaxis",
    #               label = "y-axis",
    #               value = 2,
    #               min = 1,
    #               max = 10,
    #               step = 1,
    #               round = TRUE,
    #               ticks = TRUE
    #             )
    #           )
    #         ),
    #         htmlOutput("bulk_single_cca_scores", height = "400px")
    #       ),
    #       column(
    #         width = 7,
    #         htmlOutput("bulk_single_cca", height = "500px")
    #       )
    #     )
    #   )
    # ),
  
    
    tabPanel(
      "Mass cytometry",
      
      # Show a plot of the generated distribution
      fluidPage(
        fluidRow(
          
          column(
            width = 9,
            div(
              htmlOutput("tnse_cytof", height = "400px"),
              style = "height: 400px;"
            )
          )
          
          # column(
          #   width = 3,
          #   wellPanel(radioButtons(
          #     inputId  = "cell_type",
          #     # inline   = TRUE,
          #     label    = "Cell type:",
          #     choices  = possible_cell_types_cytof,
          #     selected = "Fibroblast"
          #   ))
          # )
          
        ),
        hr(),
        fluidRow(
          column(
            width = 7,
            DT::dataTableOutput("cytof_table", height = "350px")
          )
        )
      )
    ),
    
    # tabPanel(
    #   "Clusters annotation and connection",
    #   fluidPage(
    #     h4("Table of identified single-cell RNA-seq subsets"),
    #     br(),
    #     fluidRow(
    #       column(12,
    #              dataTableOutput('table')
    #       )
    #     )
    #   )
    # ),
    
    
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
