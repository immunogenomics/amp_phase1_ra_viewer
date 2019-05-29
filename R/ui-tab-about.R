tabPanel(
  title = "About",
  value = "about",

  mainPanel(width = 12,
    
    HTML(
"
<h1>Welcome</h1>
<p>This site provides a view of some of the data from the
<a href='https://www.niams.nih.gov/grants-funding/funded-research/accelerating-medicines'>Accelerating Medicines Partnership (AMP)</a>
Rheumatoid Arthritis (RA) Phase I project.
Here, you can view single-cell RNA-seq, bulk RNA-seq, and mass cytometry data.
<p>

<div class='text-center'>
<button id='button-data' type='button' class='btn btn-primary btn-lg'>View the data</button>
</div>

<h1>Read the paper</h1>
<p>
Please cite our paper if you use this data in your work:
</p>

<div class='citation'>
<h4><a href='https://doi.org/10.1038/s41590-019-0378-1'>Defining Inflammatory Cell States in Rheumatoid Arthritis Joint Synovial Tissues by Integrating Single-cell Transcriptomics and Mass Cytometry</a></h4>
<p class='citation-authors'>Fan Zhang, Kevin Wei, Kamil Slowikowski, Chamith Y. Fonseka, Deepak A. Rao, Stephen Kelly, Susan M. Goodman, Darren Tabechian, Laura B. Hughes, Karen Salomon-Escoto, Gerald F. M. Watts, William Apruzzese, David J. Lieb, David L. Boyle, Arthur M. Mandelin II, Accelerating Medicines Partnership: RA Phase 1, AMP RA/SLE, Brendan F. Boyce, Edward DiCarlo, Ellen M. Gravallese, Peter K. Gregersen, Larry Moreland, Gary S. Firestein, Nir Hacohen, Chad Nusbaum, James A. Lederer, Harris Perlman, Costantino Pitzalis, Andrew Filer, V. Michael Holers, Vivian P. Bykerk, Laura T. Donlin, Jennifer H. Anolik, Michael B. Brenner, Soumya Raychaudhuri
</p>
<p>DOI: <a href='https://doi.org/10.1038/s41590-019-0378-1'>10.1038/s41590-019-0378-1</a></p>
<h4>Overview of synovial tissue workflow</h4>
<img src='synovial_pipeline.png' style='width:100%;'>
</div>

<h1>Download the data</h1>
<p>Download the complete data from ImmPort:</p>
<ul>
<li><a href='http://www.immport.org/immport-open/public/study/study/displayStudyDetail/SDY998'>SDY998:
AMP Rheumatoid Arthritis Arthroplasty Phase 1</a></li>
<li><a href='http://www.immport.org/immport-open/public/study/study/displayStudyDetail/SDY999'>SDY999:
AMP Rheumatoid Arthritis Synovial Phase 1</a></li>
</ul>

<h1>Get the code</h1>
<p>Get the code for data analysis:
</p>
<ul>
<li><a href='https://github.com/immunogenomics/amp_phase1_ra'>github.com/immunogenomics/amp_phase1_ra</a></li>
</ul>
<p>Get the code for this website:
</p>
<ul>
<li><a href='https://github.com/immunogenomics/amp_phase1_ra_viewer'>github.com/immunogenomics/amp_phase1_ra_viewer</a></li>
</ul>

<h1>Contact us</h1>
<p>
Please <a href='mailto:support@immunogenomics.io'>contact us</a>
if you have any questions or comments on the analysis and results.
</p>
"
    )

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
