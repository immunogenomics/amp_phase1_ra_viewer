tabPanel(
  "Home",

  mainPanel(
    h3("Welcome"),
    
    p(
      "This site provides access to view the data from the",
      a("Accelerating Medicines Partnership (AMP)",
        href = "https://www.nih.gov/research-training/accelerating-medicines-partnership-amp"),
      "Rheumatoid Arthritis (RA) Phase I project."
    ),
    br(),
    
    p(
      "To learn more about this project, please refer to our preprint paper available at",
      a("bioRxiv:",
          href = "https://www.biorxiv.org/content/early/2018/06/20/351130"),
      strong(
        "Defining Inflammatory Cell States in Rheumatoid Arthritis (RA)
        Joint Synovial Tissues by Integrating Single-cell Transcriptomics and Mass Cytometry"
      )
    ),
    br(),
    
    h5("Overview of synovial tissue workflow", align = "center"),
    img(src='synovial_pipeline.png', width = 600, align = "left"),
    br(),
    br(),

    p( 
      "Please feel free to explore the site.",
      "On this site, we provide an interactive view of the single-cell RNA-seq and mass cytometry.",
      "To access the full data, please go to",
      a("ImmPort",
        href = "http://immport.org"
      ),
      "with the study accession code of SDY998 and SDY999.",
      "The source code repository of the computational pipeline for data analysis and integration is located at",
        a("amp_phase1_ra.",
          href = "https://github.com/immunogenomics/amp_phase1_ra"
        )
    ),
    br(),
    
    p(
      "This site is maintained by", 
      a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
      "Please contact us if you have any questions, requests, or comments",
      " on the analysis and results."
    ),
    br(),
    br(),
    br()
    
    # h3("Accelerating Medicines Partnerships (AMP)"),
    # p(
    #   "The",
    #   a("Accelerating Medicines Partnership (AMP)",
    #     href = "https://www.nih.gov/research-training/accelerating-medicines-partnership-amp"),
    #   " is a public-private partnership between the National Institutes of",
    #   " Health (NIH), the U.S. Food and Drug Administration (FDA), 10",
    #   " biopharmaceutical companies and multiple non-profit organizations",
    #   " to transform the current model for developing new diagnostics and",
    #   " treatments by jointly identifying and validating promising",
    #   " biological targets for therapeutics. The ultimate goal is to",
    #   " increase the number of new diagnostics and therapies for patients",
    #   " and reduce the time and cost of developing them."
    # ),
    # 
    # h3("AMP RA Phase I"),
    # p(
    #   "Detecting distinct cellular subsets in tissues affected by rheumatoid arthritis is the key to deciphering pathogenesis in RA.",
    #   "We applied a multi-modal high dimensional strategy", 
    #   "including single-cell RNA-seq, mass cytometry, bulk RNA-seq, and flow cytometry",
    #   "to synovial tissue samples from 51 individuals with RA and osteoarthritis.",
    #   "Using an integrative strategy that uses canonical correlational analysis,", 
    #   "we are able to integrate across these data sets to define cellular populations that are robust.",
    #   "Evidence of these populations are seen across the different data modalities.",
    #   "This website supports the results of our AMP RA Phase I paper",
    #   a("(preprint version).",
    #     href = "https://www.biorxiv.org/content/early/2018/06/20/351130"),
    #     "Welcome to read the paper to know more details."
    # ),

    # h2("Disclaimer"),
    # p(
    #   "Data presented on this page is from Phase 1 of the AMP partership."
    #   # " Currently, this is private data meant to be shared internally,",
    #   # " only with consortium members."
    # ),
    # p(
    #    strong(
    #     "Sharing any data from this site with anyone outside of the",
    #    " AMP partnership is prohibited."
    #   )
    # ),
    # p(
    #   "This website is an experiment in providing early access to",
    #   " preliminary data analysis results. The content of this site is",
    #   " subject to change at any time without notice. We hope that you",
    #   " find it useful, but we provide it 'as is' without warranty of",
    #   " any kind, express or implied."
    # ),
    
    # h3("Contact"),
    # p(
    #   "This site is maintained by", 
    #   a("Kamil Slowikowski", href = "mailto:kslowikowski@fas.harvard.edu"),
    #   "and",
    #   a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
    #   "Please contact us if you have any questions, requests, or comments",
    #   " on the analysis and results."
    # ),

  ) # mainPanel
  
) # tabPanel
