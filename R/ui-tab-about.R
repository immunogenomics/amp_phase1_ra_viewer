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
      "This site is maintained by", 
      a("Kamil Slowikowski", href = "mailto:kslowikowski@fas.harvard.edu"),
      "and",
      a("Fan Zhang.", href = "mailto:fanzhang@broadinstitute.org"),
      "Please contact us if you have any questions, requests, or comments",
      " on the analysis and results."
    ),
    br()

  ) # mainPanel

) # tabPanel
