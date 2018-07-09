
message("MEMORY USAGE load-data.R 1: ", ceiling(pryr::mem_used() / 1e6), " MB")

# From BuenColors
solar_flare <- c(
  "#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", 
  "#FA8E24", "#DA2828", "#A31D1D"
)

# Differential gene expression table -------------------------------

# > head(dg)
#     gene   auc       wilcox cluster
# 1    DCN 0.733 4.858918e-57     F-1
# 2   IGF1 0.722 3.396392e-55     F-1
# 3  FBLN1 0.711 1.174966e-48     F-1
# 4 PDGFRL 0.690 1.436098e-39     F-1
# 5  SFRP1 0.688 2.526206e-41     F-1
# 6     C3 0.681 2.426171e-39     F-1
# dg_fibro <- readRDS("data/markers_gene_res_fibro.rds")
# dg_tcell <- readRDS("data/markers_gene_res_tcell.rds")
# dg_bcell <- readRDS("data/markers_gene_res_bcell.rds")
# dg_mono  <- readRDS("data/markers_gene_res_mono.rds")
# dg <- rbind(dg_fibro, dg_tcell, dg_bcell, dg_mono)
# rm(dg_fibro, dg_tcell, dg_bcell, dg_mono)
# rownames(dg) <- seq(nrow(dg))

# Upload all the genes without any filters
dg <- readRDS("data/cluster_marker_table_within_celltype_de.rds")
rownames(dg) <- seq(nrow(dg))

# > head(dg)
# gene cluster wilcox_pvalue       auc pct_nonzero
# 200229 C10orf105   SC-F4           242 0.8885902   0.8038741
# 205349     CLIC5   SC-F4           231 0.8887112   0.8111380
# 330101      SDC1   SC-B4           214 0.9617745   0.9264214
# 307561      MZB1   SC-B4           199 0.9974530   1.0000000
# 292513    FNDC3B   SC-B4           196 0.9617428   0.9397993
# 292253    FKBP11   SC-B4           195 0.9945727   0.9933110

dg$wilcox_pvalue <- round(-log10(dg$wilcox_pvalue))
dg <- dg[order(dg$wilcox_pvalue, decreasing = TRUE),]
# object_size(dg)
# 3.57 MB

# dg <- readRDS("data/cluster_marker_table.rds")
# object_size(dg)
# 78.2 MB

# Single-cell RNA-seq data -----------------------------------------

# # Read 4 datasets: bcell, tcell, mono, fibro
# # Preprocess into a file for quick loading.
# log2cpm_file <- "data/amp-phase1-ra-single-cells-matrix.h5"
# log2cpm_dimnames_file <- "data/amp-phase1-ra-single-cells-dimnames.rda"
# meta_file <- "data/amp-phase1-ra-single-cells-meta.rds"
# if (file.exists(meta_file)) {
#   meta <- readRDS(file = meta_file)
#   # ffload(log2cpm_file, overwrite = TRUE)
#   load(log2cpm_dimnames_file)
#   log2cpm <- HDF5Array::HDF5Array(file = log2cpm_file, name = "log2cpm")
#   load(log2cpm_dimnames_file)
#   rownames(log2cpm) <- log2cpm_rows
#   colnames(log2cpm) <- log2cpm_cols
# } else {
  # for (cell_type in c("fibro", "tcell", "bcell", "mono")) {
  #   assign(
  #     x = sprintf("log2cpm_%s", cell_type),
  #     value = Matrix::Matrix(
  #       data = readRDS(sprintf("data/%s_exp.rds", cell_type)),
  #       sparse = TRUE
  #     )
  #   )
  #   assign(
  #     x = sprintf("meta_%s", cell_type),
  #     value = readRDS(sprintf("data/%s_sc_label.rds", cell_type))
  #   )
  # }
  # # Combine the different cell types.
  # meta <- rbind(meta_fibro, meta_bcell, meta_tcell, meta_mono)
  # meta$cell_type <- c(
  #   rep("fibro", nrow(meta_fibro)),
  #   rep("bcell", nrow(meta_bcell)),
  #   rep("tcell", nrow(meta_tcell)),
  #   rep("mono", nrow(meta_mono))
  # )
  # meta$cell_name <- as.character(meta$cell_name)
  # meta$cluster <- as.character(meta$cluster)
  # # Change F-1 to SC-F1
  # meta$cluster <- paste("SC-", sapply(strsplit(meta$cluster, split='-', fixed=TRUE), function(x) (paste(x[1], x[2], sep=""))), sep="")
  # # Read coordinates for the all-cell tSNE plot.
  # x <- readRDS("data/all_cells_fine_cluster_label.rds")
  # # Change F-1 to SC-F1
  # x$fine_cluster <- paste("SC-", sapply(strsplit(x$fine_cluster, split='-', fixed=TRUE), function(x) (paste(x[1], x[2], sep=""))), sep="")
  # 
  # # Combine the expression matrices into one matrix.
  # log2cpm <- cbind(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
  # # Confirm that the meta data.frame matches the log2cpm matrix.
  # stopifnot(all(meta$cell_name == colnames(log2cpm)))
#   # Remove SC-T1 from log2cpm
#   log2cpm <- log2cpm[, -which(meta$cluster == "SC-T1")]
#   # Discard genes with low expression.
#   # log2cpm <- log2cpm[Matrix::rowSums(log2cpm > 0) > 10, ]
#   log2cpm <- log2cpm[rownames(log2cpm) %in% as.character(unique(dg$gene)), ]
#   
#   # Remove SC-T1 from both meta and x
#   meta <- meta[-which(meta$cluster == "SC-T1"),]
#   x <- x[-which(x$fine_cluster == "SC-T1"),]
#   
#   stopifnot(all(rownames(x) %in% meta$cell_name))
#   x <- x[meta$cell_name,]
#   meta$T1_all  <- x$T1
#   meta$T2_all  <- x$T2
#   meta$cluster <- as.character(x$fine_cluster)
#   meta$disease <- as.character(x$disease)
#   meta$plate   <- as.character(x$plate)
#   # meta$newcluster <- meta$cluster
#   # meta$newcluster[meta$cluster == "F-4"] <- "F-3"
#   # meta$newcluster[meta$cluster == "F-3"] <- "F-4"
#   # meta$cluster <- meta$newcluster
#   # meta$newcluster <- NULL
#   
#   
#   # Confirm that the meta data.frame matches the log2cpm matrix.
#   stopifnot(all(meta$cell_name == colnames(log2cpm)))
#   # Delete the temporary variables.
#   rm(log2cpm_fibro, log2cpm_bcell, log2cpm_tcell, log2cpm_mono)
#   rm(meta_fibro, meta_bcell, meta_tcell, meta_mono)
#   # Save to files.
#   # log2cpm <- as.ff(as.matrix(log2cpm), dimorder = c(2, 1))
#   # ffsave(log2cpm, file = log2cpm_file)
#   saveRDS(meta, file = meta_file)
#   log2cpm_rows <- rownames(log2cpm)
#   log2cpm_cols <- colnames(log2cpm)
#   save(list = c("log2cpm_rows", "log2cpm_cols"), file = log2cpm_dimnames_file)
#   log2cpm <- HDF5Array::writeHDF5Array(log2cpm, name = "log2cpm", file = log2cpm_file)
# }

# Use the data that we presented for the AMP RA Phase I paper
log2cpm_file <- "data/amp-phase1-ra-single-cells-matrix_5265.h5"
log2cpm_dimnames_file <- "data/amp-phase1-ra-single-cells-dimnames_5625.rda"
meta <- readRDS("data/celseq_synovium_meta_5265cells_paper.rds")
# log2cpm <- readRDS("data/celseq_synovium_log2_5265cells_paper.rds")
log2cpm <- HDF5Array::HDF5Array(file = log2cpm_file, name = "log2cpm")
stopifnot(all(meta$cell_name == colnames(log2cpm)))
load(log2cpm_dimnames_file)
rownames(log2cpm) <- log2cpm_rows
colnames(log2cpm) <- log2cpm_cols

# object_size(log2cpm)
# 2.42 MB

cluster_table <- meta %>%
  group_by(cluster) %>%
  summarise(
    cells = length(cell_name),
    RA = sum(disease == "RA"),
    OA = sum(disease == "OA")
  ) %>% as.data.frame()

gene_symbols <- rownames(log2cpm)

message("MEMORY USAGE load-data.R 2: ", ceiling(pryr::mem_used() / 1e6), " MB")

# dg_best <- dg %>% group_by(cluster) %>% top_n(n = 10, wt = -wilcox)
# auc_genes <- unique(as.character(dg_best$gene))
# mat <- log2cpm[gene_symbols %in% auc_genes,]
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
        
possible_cell_types_rna <- c(
  # "All cells"  = "all",
  # "B cell"     = "bcell",
  # "T cell"     = "tcell",
  # "Monocyte"   = "mono",
  # "Fibroblast" = "fibro"
  "All cells"  = "all",
  "B cell"     = "B cell",
  "T cell"     = "T cell",
  "Monocyte"   = "Monocyte",
  "Fibroblast" = "Fibroblast"
)

one_gene_symbol_default <- "POSTN"

meta_colors$cluster <- c(
    "SC-F1" = "#6BAED6",
    "SC-F2" = "#08306B", 
    "SC-F3" = "#DEEBF7",
    "SC-F4" = "grey",
    # "SC-T1" = "#FEB24C",
    "SC-T1" = "#8C510A",
    "SC-T2" = "brown",
    "SC-T3" = "#FFFF33",
    "SC-T4" = "#C7EAE5",
    "SC-T5" = "#003C30",
    "SC-T6" = "#35978F", 
    "SC-B1" = "#FCBBA1",
    "SC-B2" = "#CB181D", #FB6A4A #A50F15
    "SC-B3" = "#67000D",
    "SC-B4" = "#FB9A99",
    "SC-M1" = "#AE017E",
    "SC-M2" = "#F768A1",
    "SC-M3" = "#FDE0EF", #FCC5C0
    "SC-M4" = "#49006A"
)

meta_colors$cytof_cluster <- c(
  "THY1- Cadherin11-" ="#E41A1C", 
  "THY1- Cadherin11+" ="#377EB8",
  "THY1- CD34- HLA-DR+" = "#4DAF4A", 
  "THY1- CD34+ HLA-DR+" ="#984EA3", 
  "THY1+ CD34- HLA-DR-" ="#FF7F00", 
  "THY1+ CD34- HLA-DR+" ="#FFFF33",
  "THY1+ CD34+ HLA-DR-" =  "#A65628", 
  "THY1+ CD34+ HLA-DR+" ="#F781BF",
  
  "CD11c-" ="#E41A1C",   
  "CD11c+ CCR2+" ="#377EB8",
  "CD11c+ CD38-"   = "#4DAF4A",   
  "CD11c+ CD38- CD64+" ="#984EA3", 
  "CD11c+ CD38+" = "#FF7F00",
  
  "CD4- CD8-" = "#E41A1C",
  "CD4+ CCR2+" ="#377EB8",    
  "CD4+ HLA-DR+"     = "#4DAF4A", 
  "CD4+ PD-1+ ICOS-"  ="#984EA3", 
  "CD4+ PD-1+ ICOS+"  ="#FF7F00",  # 
  "CD8+ PD-1- HLA-DR-"   ="#FFFF33", # 
  "CD8+ PD-1- HLA-DR+"   =  "#A65628",    
  "CD8+ PD-1+ HLA-DR-"   ="#F781BF",
  "CD8+ PD-1+ HLA-DR+" ="#999999",
  
  "IgA+ IgM- IgD-" = "#A6CEE3",
  "IgM- IgD- HLA-DR-"="#1F78B4",
  "CD38++ CD20- IgM+ HLA-DR+" = "#B2DF8A", 
  "IgM+ IgD+ CD11c-" ="#FB9A99", 
  "IgM+ IgD+ CD11c+"  ="#E31A1C",
  "CD38++ CD20- IgM- IgD-" =  "#FDBF6F",    
  "CD38+ HLA-DR++ CD20- CD11c+" ="#FF7F00",
  "IgM+ IgD-"  ="#33A02C", 
  "IgM- IgD- HLA-DR++ CD20+ CD11c+" = "#6A3D9A",
  "IgM- IgD- HLA-DR+" = "#CAB2D6"
)

meta_colors$inflamed <- c(
  "OA" = "#6A3D9A",
  "leukocyte-poor RA" = "#FFD8B2",
  "leukocyte-rich RA" = "#FF7F00"
)

cluster_markers <- data.frame(
  Subsets = c(
    "SC-B1: IGHD+ CD27- naive B cells",
    "SC-B2: IGHG3+ CD27- memory B cells",
    "SC-B3: CD11c+ autoimmune-associated B cells",
    "SC-B4: Plasma cells",
    "SC-T1: CCR7+ CD4+ T cells",
    "SC-T2: FOXP3+ Tregs",
    "SC-T3: PD-1+ Tph/Tfh",
    "SC-T4: GZMK+ T cells",
    "SC-T5: GNLY+ GZMB+ CTLs",
    "SC-T6: GZMK+/GZMB+ T cells",
    "SC-M1: IL1B+ pro-inflammatory monocytes",
    "SC-M2: NUPR1+ monocytes",
    "SC-M3: C1QA+ monocytes",
    "SC-M4: IFN-activated monocytes",
    "SC-F1: CD34+ sublining ",
    "SC-F2: HLA+ sublining ",
    "SC-F3: DKK3+ sublining ",
    "SC-F4: CD55+ lining"
  ),
  Markers = c(
    "IGHD, CXCR4, IGHM",
    "IGHG3, HLA-DRA, MS4A1",
    "ITGAX, ACTB, TBX21, AICDA",
    "XBP1, MZB1, FKBP11, SSR4, DERL3",
    "CCR7, LEF1",
    "FOXP3, IKZF2, LAYN, CTLA4",
    "PDCD1, CXCL13",
    "GZMA, CCL5, NKG7, CD8A",
    "GNLY, GZMB, CX3CR1",
    "GZMK, HLA-DQA1, HLA-DRA",
    "IL1B, NR4A2, PLAUR, HBEGF",
    "GPNMB, HTRA1, NUPR1",
    "C1QA, MARCO",
    "SPP1, IFITM3, IFI6",
    "C3, PTGFR",
    "HLA-DRA, IL6",
    "DKK3, COL8A2",
    "CLIC5, PRG4"
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
] <- "leukocyte-rich RA"
b_meta$inflamed[
  b_meta$disease_tissue != "Arthro-OA" &
  b_meta$lymphocytes <= 0.177
] <- "leukocyte-poor RA"

b_meta$inflamed <- factor(b_meta$inflamed)
b_meta$inflamed <- fct_relevel(
  b_meta$inflamed, "OA", "leukocyte-poor RA", "leukocyte-rich RA")

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
cell_name_to_type <- meta$cell_type
names(cell_name_to_type) <- as.character(meta$cell_name)
cell_name_to_cluster <- meta$cluster
names(cell_name_to_cluster) <- as.character(meta$cell_name)
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
dat_cca$cluster <- c(
  as.character(b_meta$cell_type),
  as.character(cell_name_to_cluster[cca_bs_ynames])
)

# rbindlist(lapply(1:10, function(i) {
#   form <- as.formula(sprintf("V%s ~ cluster", i))
#   data.frame(
#     broom::tidy(anova(aov(formula = form, data = dat_cca)))[1,],
#     CV = i
#   )
# }))

# mat_cca <- data.table::dcast(
#   data = melt_cca,
#   formula = variable ~ cluster,
#   value.var = "value",
#   fun.aggregate = mean
# )
# rownames(mat_cca) <- mat_cca$variable
# mat_cca$variable <- NULL



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
# msig_names <- c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7")
# msig <- lapply(msig_names, function(i) {
#   qusage::read.gmt(sprintf("data/%s.all.v6.1.symbols.gmt", i))
# })
# names(msig) <- msig_names
# saveRDS(msig, "data/msigdb-v6.1.rds")
 
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


# Load cytof data
cytof_all <- readRDS("data/cytof_markers_tsne.rds")

# cytof_file <- "data/amp-phase1-ra-cytof-matrix.h5"
# cytof_dimnames_file <- "data/amp-phase1-ra-cytof-dimnames.rda"
# cytof_all <- HDF5Array::HDF5Array(file = cytof_file, name = "cytof_all")
# # Convert to data frame for the purpose of 
# # "cytof_all$marker <- as.numeric(cytof_all[, which(colnames(cytof_all) == marker)])" in the output
# cytof_all <- as.data.frame(cytof_all)
# load(cytof_dimnames_file)
# rownames(cytof_all) <- cytof_rows
# colnames(cytof_all) <- cytof_cols

object_size(cytof_all)
# 23.1 MB

# Take all the protein markers
protein_symbols <- colnames(cytof_all)[1:35]
colnames(cytof_all)[which(colnames(cytof_all) == "markers")] <- "cluster"

one_protein_symbol_default <- "CD90"

# Show a table of all the proteins for users to choose
proteins = data.frame(
  protein = as.character(seq(1, length(protein_symbols))),
  markers = protein_symbols
)

possible_cell_types_cytof <- c(
  "B cell"     = "B cell",
  "T cell"     = "T cell",
  "Monocyte"   = "Monocyte",
  "Fibroblast" = "Fibroblast"
)

# # Combind all the cytof cell type clusters together into one file
# 
# load("synData.Fibro.downsample.SNE.RData")
# load("synData.Bcell.downsample.SNE.RData")
# load("synData.Mono.downsample.SNE.RData")
# load("synData.Tcell.downsample.SNE.RData")
# synData.Bcell.downsample$assign_lbl <- rep("", nrow(synData.Bcell.downsample))
# 
# synData.Bcell.downsample <- synData.Bcell.downsample[, order(match(colnames(synData.Bcell.downsample), colnames(synData.Fibro.downsample)))]
# synData.Mono.downsample <- synData.Mono.downsample[, order(match(colnames(synData.Mono.downsample), colnames(synData.Fibro.downsample)))]
# synData.Tcell.downsample <- synData.Tcell.downsample[, order(match(colnames(synData.Tcell.downsample), colnames(synData.Fibro.downsample)))]
# 
# all(colnames(synData.Bcell.downsample) == colnames(synData.Fibro.downsample))
# all(colnames(synData.Bcell.downsample) == colnames(synData.Mono.downsample))
# all(colnames(synData.Bcell.downsample) == colnames(synData.Tcell.downsample))
# 
# synData.Fibro.downsample$cell_type <- rep("Fibroblast", nrow(synData.Fibro.downsample))
# synData.Mono.downsample$cell_type <- rep("Monocyte", nrow(synData.Mono.downsample))
# synData.Bcell.downsample$cell_type <- rep("B cell", nrow(synData.Bcell.downsample))
# synData.Tcell.downsample$cell_type <- rep("T cell", nrow(synData.Tcell.downsample))
# 
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90- Cadherin.11-")] <- "THY1- Cadherin11-"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90- Cadherin.11+")] <- "THY1- Cadherin11+"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90- CD34- HLA-DR+")] <- "THY1- CD34- HLA-DR+"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90- CD34+ HLA-DR+")] <- "THY1- CD34+ HLA-DR+"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90+ CD34- HLA-DR-")] <- "THY1+ CD34- HLA-DR-"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90+ CD34- HLA-DR+")] <- "THY1+ CD34- HLA-DR+"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90+ CD34+ HLA-DR-")] <- "THY1+ CD34+ HLA-DR-"
# synData.Fibro.downsample$markers[which(synData.Fibro.downsample$markers == "CD90+ CD34+ HLA-DR+")] <- "THY1+ CD34+ HLA-DR+"
# 
# synData.Tcell.downsample$markers <- as.character(synData.Tcell.downsample$markers)
# synData.Tcell.downsample$markers[which(synData.Tcell.downsample$markers == "CD8+")] <- "CD8+ PD-1- HLA-DR-"
# synData.Tcell.downsample$markers[which(synData.Tcell.downsample$markers == "CD8+ HLA-DR+")] <- "CD8+ PD-1- HLA-DR+"
# synData.Tcell.downsample$markers[which(synData.Tcell.downsample$markers == "CD8+ PD-1+")] <- "CD8+ PD-1+ HLA-DR-"
# 
# cytof_all <- rbind.data.frame(synData.Fibro.downsample, synData.Mono.downsample,
#                               synData.Tcell.downsample, synData.Bcell.downsample)
# cytof_all <- cytof_all[, c(1:35, 41:44, 46, 47, 51:53)]
# saveRDS(cytof_all, "cytof_markers_tsne.rds")
# 
# cytof_all <- as.matrix(cytof_all)
# cytof_file <- "data/amp-phase1-ra-cytof-matrix.h5"
# cytof_dimnames_file <- "data/amp-phase1-ra-cytof-dimnames.rda"
# cytof_rows <- rownames(cytof_all)
# cytof_cols <- colnames(cytof_all)
# save(list = c("cytof_rows", "cytof_cols"), file = cytof_dimnames_file)
# cytof_all <- HDF5Array::writeHDF5Array(cytof_all, name = "cytof_all", file = cytof_file)


