tabPanel(
  title = "Data Viewer",
  value = "data",
  
  h2("Clusters in single-cell RNA-seq and mass cytometry data"),
  
  fluidRow(
    column(
      width = 7,
      wellPanel(
        radioButtons(
          inputId  = "cell_type",
          inline   = TRUE,
          label    = "Cell type:",
          choices  = possible_cell_types_rna,
          selected = "all"
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
              width = 12,
              div(
                htmlOutput("tnse_marker_plot") #, height = "400px"),
                #style = "height: 400px;"
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
              DT::dataTableOutput("dg_table")#, height = "350px")
            )
          ),
        hr(),
        fluidRow(
          column(
            width = 10,
            h4("Bulk RNA-seq"),
            htmlOutput("bulk_dots", height = "500px")
          )
        )

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
            width = 12,
            div(
          selectInput(
            inputId = "cytof_marker",
            label = "Protein",
            choices = sort(c("CD45", "CD19", "RANKL", "CD64", "CD16", "CD8a", "FAP", "CD20", 
                        "CD45RO", "CD38", "PD.1", "CD14", "CD69", "CXCR5", "CD4",
                        "Podoplanin", 
                        "CD3", "CD11c", "FcRL4", "CD138", "CD90", "CCR2", "Cadherin.11", 
                        "FoxP3", "CD34", "CD146", "IgA", "TCRgd", "ICOS", "CD66b", "IgM", 
                        "VE.Cadherin", "HLA.DR", "IgD", "VCAM.1")),
            selected = "CD90",
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
              htmlOutput("tnse_cytof")
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
    )
    
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
    # br(),
    # hr(),
    
<<<<<<< Updated upstream
    mainPanel(
      h4("Cluster annotations:"),
      p(strong("SC-F1:"), "CD34+ sublining fibroblasts;",
        strong("SC-F2:"), "HLA+ sublining fibroblasts;",
        strong("SC-F3:"), "DKK3+ sublining fibroblasts;",
        strong("SC-F4:"), "CD55+ lining fibroblasts;"),
      
      p(strong("SC-M1:"), "IL1B+ pro-inflammatory monocytes;",
        strong("SC-M2:"), "NUPR1+ monocytes;",
        strong("SC-M3:"), "C1QA+ monocytes;",
        strong("SC-M4:"), "IFN-activated monocytes;"),
      
      p(strong("SC-T1:"), "CCR7+ CD4+ T cells;",
        strong("SC-T2:"), "FOXP3+ Tregs;",
        strong("SC-T3:"), "PD-1+ Tph/Tfh;",
        strong("SC-T4:"), "GZMK+ CD8+ T cells;",
        strong("SC-T5:"), "GNLY+ GZMB+ CTLs;",
        strong("SC-T6:"), "GZMK+/GZMB+ T cells;"),
      
      p(strong("SC-B1:"), "IGHD+ CD27 naive B cells;",
        strong("SC-B2:"), "IGHG3+ CD27- memory B cells;",
        strong("SC-B3:"), "autoimmune-associated cells (ABC);",
        strong("SC-B4:"), "Plasmablasts")
    ),
    br()
=======
>>>>>>> Stashed changes
    
  ), # tabsetPanel
  
    fluidRow(column(width = 12,
      h2("Cluster annotations"),
      HTML("
        <p>Fibroblasts (CD45<sup>-</sup> Podoplanin<sup>+</sup>)</p>
        <ul>
        <li>SC-F1: CD34+ sublining fibroblasts</li>
        <li>SC-F2: HLA+ sublining fibroblasts</li>
        <li>SC-F3: DKK3+ sublining fibroblasts</li>
        <li>SC-F4: CD55+ lining fibroblasts</li>
        </ul>
      
        <p>Monocytes (CD45<sup>+</sup> CD14<sup>+</sup>)</p>
        <ul>
        <li>SC-M1: IL1B+ pro-inflammatory monocytes</li>
        <li>SC-M2: NUPR1+ monocytes</li>
        <li>SC-M3: C1QA+ monocytes</li>
        <li>SC-M4: IFN-activated monocytes</li>
        </ul>
 
        <p>T cells (CD45<sup>+</sup> CD3<sup>+</sup>)</p>
        <ul>
        <li>SC-T1: CCR7+ CD4+ T cells</li>
        <li>SC-T2: FOXP3+ Tregs</li>
        <li>SC-T3: PD-1+ Tph/Tfh</li>
        <li>SC-T4: GZMK+ CD8+ T cells</li>
        <li>SC-T5: GNLY+ GZMB+ CTLs</li>
        <li>SC-T6: GZMK+/GZMB+ T cells</li>
        </ul>
 
        <p>B cells (CD45<sup>+</sup> CD3<sup>-</sup> CD19<sup>+</sup>)</p>
        <ul>
        <li>SC-B1: IGHD+ CD270 naive B cells</li>
        <li>SC-B2: IGHG3+ CD27- memory B cells</li>
        <li>SC-B3: autoimmune-associated cells (ABC)</li>
        <li>SC-B4: Plasmablasts</li>
        </ul>
      ")
    ))
  
) # tabPanel
