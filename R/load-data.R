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
    "F-1" = "#6BAED6",
    "F-2" = "#08306B",
    "F-3" = "#DEEBF7",
    "F-4" = "grey",
    "T-1" = "#FEB24C",
    "T-2" = "#8C510A",
    "T-3" = "brown",
    "T-4" = "#FFFF33",
    "T-5" = "#C7EAE5",
    "T-6" = "#003C30",
    "T-7" = "#35978F",
    "B-1" = "#FCBBA1",
    "B-2" = "#CB181D", #FB6A4A #A50F15
    "B-3" = "#67000D",
    "B-4" = "#FB9A99",
    "M-1" = "#AE017E",
    "M-2" = "#F768A1",
    "M-3" = "#FDE0EF", #FCC5C0
    "M-4" = "#49006A"
  )
)

cluster_markers <- data.frame(
  Subsets = c(
    "B-1: Naive B cells",
    "B-2: Activate B cells",
    "B-3: Age Associated B cells (ABCs)",
    "B-4: Plasma cells",
    "T-1: Naive CD4+ T cells",
    "T-2: Central memory CD4+ T cells",
    "T-3: Treg CD4+ T cells",
    "T-4: TpH/TfH CD4+ T cells",
    "T-5: GZMK+ CD8+ T cells",
    "T-6: CTL+ CD8+ T cells",
    "T-7: HLA+ CD8+ T cells",
    "M-1: activated monocytes",
    "M-2: NUPR1+",
    "M-3: C1AQ+ (common to all monocytes)",
    "M-4: IFN+",
    "F-1: THY1+ C3+ (sublining)",
    "F-2: THY1+ HLA+ (sublining)",
    "F-3: THY1+ DKK3+ (sublining)",
    "F-4: THY1- (lining)"
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

# # Change CM1 to M-1
# meta_mono$cluster <- sub("(.{1})(-*)", "\\1-\\2", substring(meta_mono$cluster, 2))
# meta_fibro$cluster <- sub("(.{1})(-*)", "\\1-\\2", substring(meta_fibro$cluster, 2))
# meta_bcell$cluster <- sub("(.{1})(-*)", "\\1-\\2", substring(meta_bcell$cluster, 2))
# meta_tcell$cluster <- sub("(.{1})(-*)", "\\1-\\2", substring(meta_tcell$cluster, 2))

