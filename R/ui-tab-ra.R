tabPanel(
  title = "Data",
  value = "data",
  
  h2("Single-cell RNA-seq, bulk RNA-seq, and mass cytometry data"),
  
  tabsetPanel(selected = "RNA-seq",
    
    tabPanel(
      "RNA-seq",
        
        # Show a plot of the generated distribution
        fluidPage(
          fluidRow(
  
            column(
              width = 12,
              div(
                radioButtons(
                  inputId  = "rnaseq_cell_type",
                  inline   = TRUE,
                  label    = "Cell type:",
                  choices  = possible_cell_types_rna,
                  selected = "all"
                ),
                htmlOutput("tnse_marker_plot") #, height = "400px"),
                #style = "height: 400px;"
              )
            )
            
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
        ),
      hr(),
      fluidRow(column(
        width = 12,
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

    ), # tabPanel

    tabPanel(
      "Mass cytometry",
      
      # Show a plot of the generated distribution
      fluidPage(
        fluidRow(
          column(
            width = 12,
            div(
              radioButtons(
                inputId  = "cytof_cell_type",
                inline   = TRUE,
                label    = "Cell type:",
                choices  = possible_cell_types_cytof,
                selected = "Fibroblast"
              ),
              selectInput(
                inputId = "cytof_marker",
                label = "Protein",
                choices = sort(c("CD45", "CD19", "RANKL", "CD64", "CD16",
                                 "CD8a", "FAP", "CD20", "CD45RO", "CD38",
                                 "PD.1", "CD14", "CD69", "CXCR5", "CD4",
                                 "Podoplanin", "CD3", "CD11c", "FcRL4",
                                 "CD138", "CD90", "CCR2", "Cadherin.11",
                                 "FoxP3", "CD34", "CD146", "IgA", "TCRgd",
                                 "ICOS", "CD66b", "IgM", "VE.Cadherin",
                                 "HLA.DR", "IgD", "VCAM.1")),
                selected = "CD90",
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL
              ),
              htmlOutput("tnse_cytof")
            )
          )
        ), # fluidRow
        hr(),
        fluidRow(
          column(
            width = 12,
            DT::dataTableOutput("cytof_table", height = "350px")
          )
        )
      ) # fluidPage
    ) # tabPanel
    
  ) # tabsetPanel
  
) # tabPanel
