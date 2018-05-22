
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
# meta$marker <- as.numeric(log2cpm[gene_ix,])
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

# fisher_file <- "data/fisher.rds"
# if (file.exists(fisher_file)) {
#   fits <- readRDS(fisher_file)
# } else {
#   fits <- Reduce(rbind, lapply(seq_along(gene_symbols), function(gene_i) {
#     marker <- gene_symbols[gene_i]
#     meta$marker <- as.numeric(log2cpm[gene_i,])
#     Reduce(rbind, lapply(clusters, function(cluster_i) {
#       fit <- broom::tidy(fisher.test(meta$cluster == cluster_i, meta$marker > 0))
#       fit$marker <- marker
#       fit$cluster <- cluster_i
#       fit
#     }))
#   }))
#   saveRDS(fits, "data/fisher.rds")
# }
# 
# # fits %>% arrange(p.value) %>% head()
# 
# fits$p.value[fits$p.value == 0] <- 1
# 
# fits %>% group_by(cluster) %>% top_n(n = 3, wt = -log10(p.value))

# clusters <- sort(unique(as.character(meta$cluster)))

# library(parallel)
# library(lme4)
# 
# clusters <- sprintf("M-%s", 1:4)
# ix_m <- which(meta$cluster %in% clusters)
# meta_m <- meta[ix_m,]
# 
# fits <- Reduce(rbind, mclapply(seq_along(gene_symbols)[1:100], function(gene_i) {
#   marker <- gene_symbols[gene_i]
#   meta_m$marker <- as.numeric(log2cpm[gene_i,ix_m])
#   Reduce(rbind, lapply(clusters, function(cluster_i) {
#     # fit <- broom::tidy(fisher.test(meta$cluster == cluster_i, meta$marker > 0))
#     fit <- broom::tidy(glmer(cluster == cluster_i ~ marker + (1|plate), meta_m, family = binomial))
#     fit$marker <- marker
#     fit$cluster <- cluster_i
#     fit
#   }))
# }, mc.cores = 4))
# 
# fits$p.value[fits$p.value == 0] <- 1
# fits$p.value[is.na(fits$p.value)] <- 1
# 
# x <- fits[fits$term == "marker",]
# x %>% group_by(cluster) %>% top_n(n = 5, wt = -log10(p.value))

# gene_i <- which(gene_symbols == "ABI3")
# meta$marker <- as.numeric(log2cpm[gene_i,])
# auroc1l(meta$marker, meta$cluster == "M-2")

# # https://gist.github.com/mbq/21ef2370961a634ce45b44007e938b60
# auroc <- function(score, cls) {
#   n1 <- sum(!cls)
#   n2 <- sum(cls)
#   U <- sum(rank(score)[!cls]) - n1 * (n1 + 1) / 2
#   return(1 - U / n1 / n2)
# }
# # p-value of having AUROC of auroc or more
# # with nx positive and ny negative classes
# auroc_p <-function(auroc, nx, ny) {
#   W <- round((1 - auroc) * nx * ny);
#   pwilcox(W, nx, ny)
# }

# 
# 
# auc_m <- Reduce(rbind, lapply(gene_symbols_m, function(gene_symbol) {
#   gene_i <- which(gene_symbols == gene_symbol)
#   meta_m$marker <- as.numeric(log2cpm[gene_i, ix_m])
#   sapply(X = clusters, FUN = function(cluster_i) {
#     indicator <- meta_m$cluster == cluster_i
#     auc <- auroc(meta_m$marker, indicator)
#     auc_p <- auroc_p(auc, sum(indicator), sum(!indicator))
#     c(auc, auc_p)
#   })
# }))
# rownames(auc_m) <- gene_symbols_m
# 
# auc_m[auc_m < 0.5] <- 1 - auc_m[auc_m < 0.5]
# 
# apply(auc_m, 2, function(col) {
#   which(col > 0.77)
# })
# 
# diff_m <- apply(auc_m, 1, function(row) {
#   xs <- sort(row, decreasing = TRUE)
#   xs[1] - xs[2]
# })


library(MUDAN)
library(reshape2)
library(pheatmap)
library(BuenColors)

# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

mudan_file <- "data/mudan_differential_genes.rda"
if (file.exists(mudan_file)) {
  load(mudan_file)
} else {
  get_markers <- function(clusters) {
    # clusters <- sprintf("M-%s", 1:4)
    ix_m <- which(meta$cluster %in% clusters)
    meta_m <- meta[ix_m,]
    genes_m <- vapply(seq_along(gene_symbols), function(gene_i) {
      x <- as.numeric(log2cpm[gene_i,ix_m])
      sum(x > 0)
    }, 1.0)
    gene_symbols_m <- gene_symbols[genes_m > 10]
    ix_genes <- which(gene_symbols %in% gene_symbols_m)
    mat_m <- as.numeric(log2cpm[ix_genes, ix_m])
    rownames(mat_m) <- gene_symbols[ix_genes]
    colnames(mat_m) <- meta_m$cell_name
    cols <- as.character(meta_m$cluster)
    names(cols) <- meta_m$cell_name
    d <- MUDAN::getDifferentialGenes(mat_m, cols)
    for (i in 1:length(d)) {
      d[[i]]["cluster"] <- names(d)[i]
      d[[i]]["gene"] <- rownames(d[[i]])
      rownames(d[[i]]) <- seq(nrow(d[[i]]))
    }
    d <- do.call(rbind, d)
    rownames(d) <- seq(nrow(d))
    return(d)
  }
  
  d_m <- get_markers(sprintf("M-%s", 1:4))
  d_f <- get_markers(sprintf("F-%s", 1:4))
  d_t <- get_markers(sprintf("T-%s", 1:7))
  d_b <- get_markers(sprintf("B-%s", 1:4))
  
  save(
    list = c("d_m", "d_f", "d_t", "d_b"),
    file = mudan_file
  )
}

plot_z <- function(d) {
  ggplot() +
    geom_point(
      data = d[d$highest,],
      mapping = aes(x = fe * 100, y = Z),
      size = 0.5
    ) +
    labs(x = "Percent of nonzero cells") +
    facet_grid(~ cluster) +
    theme_bw(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(1, "lines")
    )
}

plot_z(d_m)
ggsave("d_m.png", width = 8, height = 3)
plot_z(d_f)
ggsave("d_f.png", width = 8, height = 3)
plot_z(d_t)
ggsave("d_t.png", width = 12, height = 3)
plot_z(d_b)
ggsave("d_b.png", width = 8, height = 3)

plot_heat <- function(d, n = 5) {
  best_genes <- subset(d, highest == TRUE & fe > 0.75) %>%
    group_by(cluster) %>%
    top_n(n = n, wt = Z)
  
  mat <- dcast(
    data = d[d$gene %in% best_genes$gene,],
    formula = gene ~ cluster,
    value.var = "Z"
  )
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  
  pheatmap(
    mat = mat,
    color = colorRampPalette(BuenColors::jdb_palettes$solar_flare)(20),
    border_color = NA, fontsize = 20
  )
}

png("d_m_heat.png", width = 400, height = 600)
plot_heat(d_m)
dev.off()

png("d_f_heat.png", width = 400, height = 600)
plot_heat(d_f)
dev.off()

png("d_t_heat.png", width = 400, height = 600)
plot_heat(d_t)
dev.off()

png("d_b_heat.png", width = 400, height = 600)
plot_heat(d_b)
dev.off()
