# Differential gene expression table.
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

# 
# library(lme4) # needed for `VarCorr` function
# library(lme4qtl)
# 
# # load synthetic data set `dat40` distributed within `lme4qtl`
# # - table of phenotypes `dat40`
# # - the double kinship matrix `kin2`
# data(dat40)
# 
# # (1) model continiuous trait `trait1`
# mod <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))
# 
# # get the estimation of h2
# (vf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")])
# #       grp      vcov
# #1       ID 5.2845001
# #2    FAMID 0.0000000
# #3 Residual 0.6172059
# 
# prop <- with(vf, vcov / sum(vcov))
# 
# (h2 <- prop[1]) 
# #[1] `0.895419`
# 
# # (2) model binary trait `trait1bin`
# gmod <- relmatGlmer(trait1bin ~ (1|ID), dat40, relmat = list(ID = kin2), family = binomial)
# 
# marker <- "IL6"
# gene_ix <- which(gene_symbols == marker)
# meta$marker <- lf$matrix[,gene_ix]
# 
# clusters <- sort(unique(as.character(meta$cluster)))
# fits <- Reduce(rbind, mclapply(clusters, function(cluster_i) {
#   fit <- broom::tidy(glmer(cluster == cluster_i ~ marker + (1|plate), meta, family = binomial))
#   fit$marker <- marker
#   fit$cluster <- cluster_i
#   fit
# }, mc.cores = 4))
# fits$p.value[fits$p.value == 0] <- 1
# 
# fits <- Reduce(rbind, lapply(seq_along(gene_symbols), function(gene_i) {
#   marker <- gene_symbols[gene_i]
#   meta$marker <- lf$matrix[,gene_i]
#   Reduce(rbind, lapply(clusters, function(cluster_i) {
#     fit <- broom::tidy(fisher.test(meta$cluster == cluster_i, meta$marker > 0))
#     fit$marker <- marker
#     fit$cluster <- cluster_i
#     fit
#   }))
# }))
# saveRDS(fits, "data/fisher.rds")
# 
# fits %>% arrange(p.value) %>% head()
# 
# fits$p.value[fits$p.value == 0] <- 1



