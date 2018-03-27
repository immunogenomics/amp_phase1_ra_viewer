
# Differential gene expression table -------------------------------

# > head(dg)
#     gene   auc       wilcox cluster
# 1    DCN 0.733 4.858918e-57     F-1
# 2   IGF1 0.722 3.396392e-55     F-1
# 3  FBLN1 0.711 1.174966e-48     F-1
# 4 PDGFRL 0.690 1.436098e-39     F-1
# 5  SFRP1 0.688 2.526206e-41     F-1
# 6     C3 0.681 2.426171e-39     F-1
dg_fibro <- readRDS("data/markers_gene_res_fibro.rds")
dg_tcell <- readRDS("data/markers_gene_res_tcell.rds")
dg_bcell <- readRDS("data/markers_gene_res_bcell.rds")
dg_mono  <- readRDS("data/markers_gene_res_mono.rds")
dg <- rbind(dg_fibro, dg_tcell, dg_bcell, dg_mono)
rm(dg_fibro, dg_tcell, dg_bcell, dg_mono)
rownames(dg) <- seq(nrow(dg))
# dg$newcluster <- dg$cluster
# dg$newcluster[dg$cluster == "F-4"] <- "F-3"
# dg$newcluster[dg$cluster == "F-3"] <- "F-4"
# dg$cluster <- dg$newcluster
# dg$newcluster <- NULL
dg$wilcox <- round(-log10(dg$wilcox))
dg <- dg[order(dg$wilcox, decreasing = TRUE),]

# Single-cell RNA-seq data -----------------------------------------

# Read 4 datasets: bcell, tcell, mono, fibro
# Preprocess into a file for quick loading.
loom_file <- "data/amp-phase1-ra-single-cells.loom"
if (file.exists(loom_file)) {
  lf <- loomR::connect(filename = loom_file, mode = "r", skip.validate = FALSE)
} else {
  for (cell_type in c("fibro", "tcell", "bcell", "mono")) {
    assign(
      x = sprintf("log2cpm_%s", cell_type),
      value = Matrix::Matrix(
        data = readRDS(sprintf("data/%s_exp.rds", cell_type)),
        sparse = TRUE
      )
    )
    assign(
      x = sprintf("meta_%s", cell_type),
      value = readRDS(sprintf("data/%s_sc_label.rds", cell_type))
    )
  }
  # Combine the different cell types.
  meta <- rbind(meta_fibro, meta_bcell, meta_tcell, meta_mono)
  meta$cell_type <- c(
    rep("fibro", nrow(meta_fibro)),
    rep("bcell", nrow(meta_bcell)),
    rep("tcell", nrow(meta_tcell)),
    rep("mono", nrow(meta_mono))
  )
  meta$cell_name <- as.character(meta$cell_name)
  meta$cluster <- as.character(meta$cluster)
  # Read coordinates for the all-cell tSNE plot.
  x <- readRDS("data/all_cells_fine_cluster_label.rds")
  stopifnot(all(rownames(x) %in% meta$cell_name))
  x <- x[meta$cell_name,]
  meta$T1_all  <- x$T1
  meta$T2_all  <- x$T2
  meta$cluster <- as.character(x$fine_cluster)
  meta$disease <- as.character(x$disease)
  meta$plate   <- as.character(x$plate)
  # meta$newcluster <- meta$cluster
  # meta$newcluster[meta$cluster == "F-4"] <- "F-3"
  # meta$newcluster[meta$cluster == "F-3"] <- "F-4"
  # meta$cluster <- meta$newcluster
  # meta$newcluster <- NULL
  # Combine the expression matrices into one matrix.
  log2cpm <- cbind(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
  # Discard genes with low expression.
  # log2cpm <- log2cpm[Matrix::rowSums(log2cpm > 0) > 10, ]
  log2cpm <- log2cpm[rownames(log2cpm) %in% as.character(unique(dg$gene)), ]
  # Confirm that the meta data.frame matches the log2cpm matrix.
  stopifnot(all(meta$cell_name == colnames(log2cpm)))
  # Delete the temporary variables.
  rm(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
  rm(meta_fibro, meta_bcell, meta_tcell, meta_mono)
  # Create a loom file for quick and easy gene lookups in the app.
  lf <- loomR::create(
    filename   = loom_file,
    data       = Matrix::t(log2cpm), # rows are cells, columns are genes
    cell.attrs = meta
  )
}

meta <- lf$col.attrs
col_names <- names(meta)
meta <- as.data.frame(lapply(col_names, function(col_name) {
  meta[[col_name]][]
}))
colnames(meta) <- col_names

cluster_table <- meta %>%
  group_by(cluster) %>%
  summarise(
    cells = length(cell_name),
    RA = sum(disease == "RA"),
    OA = sum(disease == "OA")
  ) %>% as.data.frame()

gene_symbols <- lf$row.attrs$gene_names[]

# dg_best <- dg %>% group_by(cluster) %>% top_n(n = 10, wt = -wilcox)
# auc_genes <- unique(as.character(dg_best$gene))
# mat <- lf$matrix[,gene_symbols %in% auc_genes]
# meta2 <- cbind(meta[,c("cluster"),drop=FALSE], mat)
# 
# cluster_mat <- reshape2::melt(meta2, id.vars = "cluster") %>%
#   group_by(cluster, variable) %>% summarise(median = median(value))
# cluster_mat <- data.table::dcast(cluster_mat, cluster ~ variable, value.var = "median")
# rownames(cluster_mat) <- cluster_mat$cluster
# cluster_mat$cluster <- NULL
# colnames(cluster_mat) <- auc_genes
# cluster_mat <- as.matrix(cluster_mat)
# dim(cluster_mat)
# 
# # cluster_mat_ix <- (colMeans(cluster_mat) > 3) &
# #   (colSums(cluster_mat > 0) < nrow(cluster_mat) / 2)
# # cluster_mat <- cluster_mat[,colMeans(cluster_mat) > 3] > 0
# # cluster_mat <- cluster_mat[,colSums(cluster_mat) < nrow(cluster_mat) / 2]
# # hmap(cluster_mat[,cluster_mat_ix], method = "TSP")
# 
# cluster_mat_ix <- (colMeans(cluster_mat) > 1)
# mat <- t(cluster_mat[,cluster_mat_ix])
# # mat <- t(cluster_mat)
# x <- seriate(mat, method = "BEA_TSP")
# 
# pheatmap::pheatmap(
#   mat = mat[x[[1]], x[[2]]],
#   cluster_cols = FALSE,
#   cluster_rows = FALSE,
#   border_color = NA
# )
        
cell_types   <- lf$col.attrs$cell_type[]

possible_cell_types <- c(
  "All cells"  = "all",
  "B cell"     = "bcell",
  "T cell"     = "tcell",
  "Monocyte"   = "mono",
  "Fibroblast" = "fibro"
)

# one_gene_symbol_default <- "HLA-DRA"
one_gene_symbol_default <- "POSTN"

meta_colors$fine_cluster <- c(
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
meta_colors$inflamed <- c(
  "OA" = "#6A3D9A",
  "RA" = "#FFD8B2",
  "inflamed RA" = "#FF7F00"
)

cluster_markers <- data.frame(
  Subsets = c(
    "B-1: Naive B cells",
    "B-2: Memory B cells",
    "B-3: Age-associated B cells/activated B cells",
    "B-4: Plasma cells",
    "T-1: Effector memory T cells",
    "T-2: Central memory T cells",
    "T-3: Treg",
    "T-4: Tph/Tfh",
    "T-5: GZMK+ CD8+ T cells",
    "T-6: CTL+ CD8+ T cells",
    "T-7: HLA+ CD8+ T cells",
    "M-1: pro-inflammatory",
    "M-2: NUPR1+",
    "M-3: C1QA+ (common to all monocytes)",
    "M-4: IFN-activated",
    "F-1: Sublining CD34+",
    "F-2: Sublining HLA+ IFN+",
    "F-3: Sublining DKK3+",
    "F-4: Lining"
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
    "CLIC5, PRG4",
    "DKK3, COL8A2"
  )
)

# Bulk RNA-seq data ------------------------------------------------

b_log2tpm <- readRDS("data/filtered_log2tpm_lowinput_phase_1.rds")
b_log2tpm <- as.matrix(b_log2tpm)
b_meta <- readRDS("data/filtered_meta_lowinput_phase_1.rds")
b_meta <- janitor::clean_names(b_meta)
stopifnot(all(b_meta$Sample.ID == colnames(b_log2tpm)))

b_meta$cell_type <- factor(
  b_meta$cell_type, rev(c("Mono", "Fibro", "B cell", "T cell"))
)

# Fix typos
b_meta$cell_type[b_meta$sample_id == "S81"] <- "Mono"
b_meta$cell_type[b_meta$sample_id == "S82"] <- "Fibro"
b_meta$cell_type[b_meta$sample_id == "S163"] <- "B cell"
b_meta$cell_type[b_meta$sample_id == "S164"] <- "T cell"

b_meta$inflamed <- "OA"
b_meta$inflamed[
  b_meta$disease_tissue != "Arthro-OA" &
  b_meta$lymphocytes > 0.177
] <- "inflamed RA"
b_meta$inflamed[
  b_meta$disease_tissue != "Arthro-OA" &
  b_meta$lymphocytes <= 0.177
] <- "RA"

b_meta$inflamed <- factor(b_meta$inflamed)
b_meta$inflamed <- fct_relevel(
  b_meta$inflamed, "inflamed RA", "RA", "OA")

# b_genes <- data.frame(
#   mean = rowMeans(b_log2tpm),
#   sd = rowSds(b_log2tpm)
# )
# ggplot() +
#   geom_point(
#     data = b_genes,
#     mapping = aes(x = mean, y = sd),
#     size = 0.5
#   ) +
#   theme_clean()
# b_mat <- scale_rows(scale(b_log2tpm))

# table(b_meta$inflamed)

# marker <- "CLIC5"
# b_meta$marker <- as.numeric(b_log2tpm[marker,])
# plot_bulk_dots(b_meta, "CLIC5")

# d <- subset(b_meta, auto_calculation_of_das28_crp > 0)
# numeric_cols <- sapply(colnames(d), function(x) {
#   is.numeric(d[[x]])
# })
# numeric_cols <- names(numeric_cols)[numeric_cols]
# fits <- lapply(numeric_cols, function(x) {
#   wilcox.test(d[["auto_calculation_of_das28_crp"]], d[[x]])
# })

# CCA between bulk and single-cell RNA-seq -------------------------

cca_bs <- readRDS("data/allcelltypes_cca_7465_genes.rds")

# 167 pairs of canonical variates

stopifnot(all(cca_bs$names$Xnames == b_meta$sample_id))
cca_bs_xnames <- cca_bs$names$Xnames

# TODO The CCA analysis should use the same cells
# as differential expression analysis.
# all(cca_bs$names$Ynames == meta$cell_name)

ix <- which(cca_bs$names$Ynames %in% meta$cell_name)
all(cca_bs$names$Ynames[ix] %in% meta$cell_name)
cca_bs_ynames <- cca_bs$names$Ynames[ix]

dat_cca <- as.data.frame(rbind(
  cca_bs$scores$corr.X.xscores[,1:10],
  cca_bs$scores$corr.Y.yscores[,1:10]
))
dat_cca <- dat_cca[c(cca_bs_xnames, cca_bs_ynames),]
cell_name_to_type <- structure(
  .Data = meta$cell_type,
  .Names = as.character(meta$cell_name)
)
cell_name_to_type <- fct_recode(
  cell_name_to_type,
  "Fibro" = "fibro",
  "B cell" = "bcell",
  "T cell" = "tcell",
  "Mono" = "mono"
)
dat_cca$cell_type <- c(
  as.character(b_meta$cell_type),
  as.character(cell_name_to_type[cca_bs_ynames])
)
# table(dat_cca$cell_type)
dat_cca$data <- c(
  rep("bulk", nrow(b_meta)),
  rep("cell", nrow(dat_cca) - nrow(b_meta))
)

# Discover mislabeled bulk samples
# ggplot() +
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_vline(xintercept = 0, color = "grey90") +
#   geom_point(
#     data = subset(dat_cca, data == "bulk"),
#     mapping = aes(
#       x = V3, y = V4, fill = cell_type,
#       shape = data,
#       size = data
#     ),
#     stroke = 0.1
#   ) +
#   scale_shape_manual(values = c(22, 21)) +
#   scale_size_manual(values = c(3, 2)) +
#   # geom_circle(
#   #   mapping = aes(x0 = 0, y0 = 0, r = 1)
#   # ) +
#   coord_equal() +
#   facet_wrap(~ cell_type) +
#   theme_clean(base_size = 20)
# subset(
#   dat_cca,
#   data == "bulk" & V3 < -0.2 & cell_type %in% c("B cell", "Mono")
# )
# subset(
#   dat_cca,
#   data == "bulk" & V3 > -0.1 & !cell_type %in% c("B cell", "Mono")
# )
# weird <- c("S81", "S82", "S163", "S164")

# Gene set enrichment ----------------------------------------------

# library(qusage)
# msigdb_hallmark <- qusage::read.gmt("data/h.all.v6.1.symbols.gmt")
# 
# library(liger)
# lig <- liger::bulk.gsea(
#   values = cca_bs$scores$xscores[,1],
#   set.list = msigdb_hallmark
# )

# > str(cca_bs)
# List of 5
# $ cor   : num [1:167] 1 1 1 1 1 ...
# $ names :List of 3
# ..$ Xnames   : chr [1:167] "S232" "S10" "S18" "S26" ...
# ..$ Ynames   : chr [1:7127] "S006_L1Q1_A01" "S006_L1Q1_A03" "S006_L1Q1_A05" "S006_L1Q1_A07" ...
# ..$ ind.names: chr [1:7465] "DPM1" "CFH" "FUCA2" "ANKIB1" ...

# xcoef has the coefficients on each bulk sample
# ycoef has the coefficients on each single cell

# $ xcoef : num [1:167, 1:167] -0.3077 0.0333 0.1558 -0.0738 -0.3118 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:167] "S232" "S10" "S18" "S26" ...
# .. ..$ : NULL
# $ ycoef : num [1:7127, 1:167] -0.0256 0.0271 -0.2266 -0.189 0.1305 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:7127] "S006_L1Q1_A01" "S006_L1Q1_A03" "S006_L1Q1_A05" "S006_L1Q1_A07" ...
# .. ..$ : NULL

# xscores has the positions of bulk genes along the 167 CVs

# $ scores:List of 6
# ..$ xscores       : num [1:7465, 1:167] 0.688 0.072 0.314 -0.861 -0.229 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:7465] "DPM1" "CFH" "FUCA2" "ANKIB1" ...
# .. .. ..$ : NULL

# yscores has the positions of single-cell genes along the 167 CVs

# ..$ yscores       : num [1:7465, 1:167] 0.6915 0.0715 0.3147 -0.8599 -0.2323 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:7465] "DPM1" "CFH" "FUCA2" "ANKIB1" ...
# .. .. ..$ : NULL

# rows of `corr.X.xscores` have the correlations of the bulk
# samples with the canonical variates `xscores`

# ..$ corr.X.xscores: num [1:167, 1:167] 0.733 0.616 0.593 0.573 0.603 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:167] "S232" "S10" "S18" "S26" ...
# .. .. ..$ : NULL
# ..$ corr.Y.xscores: num [1:7127, 1:167] 0.215 0.328 0.363 0.273 0.31 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:7127] "S006_L1Q1_A01" "S006_L1Q1_A03" "S006_L1Q1_A05" "S006_L1Q1_A07" ...
# .. .. ..$ : NULL
# ..$ corr.X.yscores: num [1:167, 1:167] 0.733 0.616 0.593 0.573 0.603 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:167] "S232" "S10" "S18" "S26" ...
# .. .. ..$ : NULL
# ..$ corr.Y.yscores: num [1:7127, 1:167] 0.215 0.328 0.363 0.273 0.31 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:7127] "S006_L1Q1_A01" "S006_L1Q1_A03" "S006_L1Q1_A05" "S006_L1Q1_A07" ...
# .. .. ..$ : NULL