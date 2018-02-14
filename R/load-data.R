# Read 4 datasets: bcell, tcell, mono, fibro
# Preprocess into a file for quick loading.
data_file <- "data/shiny.rda"
# if (file.exists(data_file)) {
#   load(data_file)
# } else {
#   e <- environment()
#   for (cell_type in cell_types) {
#     log2cpm <- readRDS(sprintf("data/%s_exp.rds", cell_type))
#     log2cpm <- Matrix(log2cpm)
#     meta    <- readRDS(sprintf("data/%s_sc_label.rds", cell_type))
#     nonzero <- rownames(log2cpm)[
#        which(rowSums(log2cpm > 0) > 10)
#     ]
#     # log2cpm_filter <- log2cpm[nonzero,]
#     # markers <- get_markers(log2cpm_filter, meta$cluster)
#     assign(sprintf("log2cpm_%s", cell_type), log2cpm, envir = e)
#     assign(sprintf("meta_%s", cell_type), meta, envir = e)
#     assign(sprintf("nonzero_%s", cell_type), nonzero, envir = e)
#     # assign(sprintf("markers_%s", cell_type), markers, envir = e)
#   }
#   #all_cell_types <- cbind.data.frame(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
#   all_meta <- rbind.data.frame(meta_fibro, meta_bcell, meta_tcell, meta_mono)
#   rm(log2cpm)
#   rm(meta)
#   rm(nonzero)
#   # rm(markers)
#   gene_symbols <- unique(c(
#     nonzero_bcell, nonzero_fibro, nonzero_mono, nonzero_tcell
#   ))
#   save.image(data_file)
# }

loom_file <- "data/amp-phase1-ra-single-cells.loom"
if (file.exists(loom_file)) {
  lf <- loomR::connect(filename = loom_file, mode = "r", skip.validate = FALSE)
  meta <- lf$col.attrs
  col_names <- names(meta)
  meta <- as.data.frame(lapply(col_names, function(col_name) {
    meta[[col_name]][]
  }))
  colnames(meta) <- col_names
} else {
  load(data_file)
  log2cpm <- cbind(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
  log2cpm <- log2cpm[rowSums(log2cpm > 0 ) > 10,]
  n_cells <- ncol(log2cpm)
  n_genes <- nrow(log2cpm)
  meta    <- rbind.data.frame(meta_fibro, meta_bcell, meta_tcell, meta_mono)
  stopifnot(all(colnames(log2cpm) == rownames(meta)))
  meta$cell_type <- c(
    rep("fibro", nrow(meta_fibro)),
    rep("bcell", nrow(meta_bcell)),
    rep("tcell", nrow(meta_tcell)),
    rep("mono", nrow(meta_mono))
  )
  meta$cell_name <- as.character(meta$cell_name)
  meta$cluster <- as.character(meta$cluster)
  lf <- loomR::create(
    filename   = loom_file,
    data       = t(log2cpm), # rows are cells, columns are genes
    cell.attrs = meta
  )
}

gene_symbols <- lf$row.attrs$gene_names[]

cell_types   <- lf$col.attrs$cell_type[]

possible_cell_types <- c(
  "B cell"     = "bcell",
  "T cell"     = "tcell",
  "Monocyte"   = "mono",
  "Fibroblast" = "fibro"
)

one_gene_symbol_default <- "HLA-DRA"

meta_colors <- list(
  fine_cluster = c(
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

cluster_markers <- data.frame(
  Subsets = c(
    "CB1 (Naive B cells)",
    "CB2 (Activate B cells)",
    "CB3 (ABCs)",
    "CB4 (Plasma cells)",
    "CT1 (Naive CD4+ T cells)",
    "CT2 (Central memory CD4+ T cells)",
    "CT3 (Treg CD4+ T cells)",
    "CT4 (TpH/TfH CD4+ T cells)",
    "CT5 (GZMK CD8+ T cells)",
    "CT6 (CTL CD8+ T cells)",
    "CT7 (HLA CD8+ T cells)",
    "CM1 (PLAUR+)",
    "CM2 (NUPR1+)",
    "CM3 (C1AQ+)",
    "CM4 (IFN+)",
    "CF1 (THY1+ C3+)",
    "CF2 (THY1+ HLA+)",
    "CF3 (THY1+ DKK3+)",
    "CF4 (THY1-)"
  ),
  Markers = c(
    "IGHD, CXCR4, IGHM",
    "HLA-DPB1, HLA-DRA, MS4A1",
    "ITGAX, ACTB, TBX21, AICDA",
    "XBP1, MZB1, FKBP11, SSR4, DERL3",
    "PTPRC",
    "CCR7, LEF1",
    "FOXP3, IKZF2, LAYN, CTLA4",
    "PDCD1, CXCL13",
    "GZMA, CCL5, NKG7, CD8A",
    "GNLY, CX3CR1, GZMB",
    "HLA-DQA1, HLA-DRA",
    "NR4A2, PLAUR, HBEGF",
    "GPNMB, HTRA1, NUPR1",
    "C1QA, MARCO",
    "SPP1, IFITM3, IFI6",
    "C3, PTGFR",
    "HLA-DRA, IL6",
    "DKK3, COL8A2",
    "CLIC5, PRG4"
  )
)
