# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'writexl', 'here', 'fpc','NbClust')
for (i in libraries) {
  library(i, character.only = TRUE)
}

# directory definition *********************************************************
cbio_dir <- Sys.getenv("CBIO_DATA")
tcga_dir <- Sys.getenv("TCGA_DATA")
pub_stad_dir <- file.path(cbio_dir, "stad_tcga_pub")
gdc_stad_dir <- file.path(cbio_dir, "stad_tcga_gdc")


# Retreiving cBioPortal data 1**************************************************
cbiopub_clin_sample <- read.delim(file.path(pub_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
cbiopub_clin_pat <- read.delim(file.path(pub_stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)
# Retreiving cBioPortal data 2**************************************************
cbiogdc_clin_sample <- read.delim(file.path(gdc_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
# Retreiving TCGA data from .rds files *****************************************
clinical_tcga <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_clinical.rds"))
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))


# ******************************************************************************
# CLUSTERING
# ******************************************************************************
# data retrieving from 01_mRNA_data.R
tpm_ind <- readRDS(here("misc", "data", "tpm_ind.rds"))
# t(matrix) para que la disimilitud se calcule usando los genes como variables
tpm_z <- t(scale(t(tpm_ind)))


# ******************************************************************************
# clustering 1: euclidean + ward.D2
# ******************************************************************************
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")

# ******************************************************************************
# clustering 2: 1-Pearson + average
# ******************************************************************************
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
hcl_p <- hclust(dist_pearson, method = "average")

# ******************************************************************************
# corte - asociación de muestras a clústers
# ******************************************************************************

k <- 3

cluster_hcl <- factor(cutree(hcl, k), levels = 1:3)
cluster_hcl_p <- factor(cutree(hcl_p, k), levels = 1:3)

# ******************************************************************************
# genes diferencialmente expresados entre clústers HCL
# ******************************************************************************
counts <- assay(rnaseq_se, "unstranded")
counts <- counts[, names(cluster_hcl), drop = FALSE]
stopifnot(identical(colnames(counts), names(cluster_hcl)))

cluster <- factor(cluster_hcl[colnames(counts)], levels = 1:k, labels = paste0("C", 1:k))
dge <- DGEList(counts = counts, group = cluster)
keep <- filterByExpr(dge, group = cluster)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + cluster)
colnames(design) <- levels(cluster)

v <- voom(dge, design)
fit <- lmFit(v, design)

con <- paste0("C", seq_len(k))
mc <- character(k)

for (i in seq_along(con)) {
  mc[i] <- paste0(con[i], " - (", paste(con[-i], collapse = " + "), ") / ", k - 1)
}

contrasts_mat <- makeContrasts(contrasts = mc, levels = design)
colnames(contrasts_mat) <- paste0(con, "all")

fit_clusters <- contrasts.fit(fit, contrasts_mat)
fit_ebayes <- eBayes(fit_clusters)

genes_C1 <- topTable(fit_ebayes, coef = "C1all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C2 <- topTable(fit_ebayes, coef = "C2all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C3 <- topTable(fit_ebayes, coef = "C3all", number = Inf, p.value = 0.05, sort.by = "logFC")

info_C1 <- head(genes_C1, 10)
info_C2 <- head(genes_C2, 10)
info_C3 <- head(genes_C3, 10)


# ******************************************************************************
# heatmap
# ******************************************************************************
gene_info <- as.data.frame(rowData(rnaseq_se))

info_C1$gene_id <- rownames(info_C1)
info_C2$gene_id <- rownames(info_C2)
info_C3$gene_id <- rownames(info_C3)

info_C1$gene_symbol <- gene_info$gene_name[match(info_C1$gene_id, rownames(gene_info))]
info_C2$gene_symbol <- gene_info$gene_name[match(info_C2$gene_id, rownames(gene_info))]
info_C3$gene_symbol <- gene_info$gene_name[match(info_C3$gene_id, rownames(gene_info))]

info_C3$gene_cluster <- "C1"
info_C2$gene_cluster <- "C2"
info_C1$gene_cluster <- "C3"

columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")
info_C1 <- info_C1[, columns_keep]
info_C2 <- info_C2[, columns_keep]
info_C3 <- info_C3[, columns_keep]

genes_heatmap_info <- rbind(info_C1, info_C2, info_C3)

gene_ids <- genes_heatmap_info$gene_id

mat_heatmap <- v$E[gene_ids, , drop = FALSE]

mat_z <- mat_heatmap
for (i in 1:nrow(mat_heatmap)) {
  gene_mean <- mean(mat_heatmap[i, ])
  gene_sd <- sd(mat_heatmap[i, ])
  mat_z[i, ] <- (mat_heatmap[i, ] - gene_mean)/gene_sd
}

gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = c("C1", "C2", "C3"))
sample_cluster <- factor(cluster[colnames(mat_z)], levels = c("C1", "C2", "C3"))

anyDuplicated(genes_heatmap_info$gene_id)

cluster_colors <- c("C1" = "#FF0000","C2" = "blue","C3" = "#00FF00")

annotation_clusters <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors),
    labels = c("C1", "C2", "C3"),
    labels_gp = gpar(col = "white", fontface = "bold")),
  which = "column")

sample_id_mat <- substr(colnames(mat_z), 1, 15)
posicion <- match(sample_id_mat, cbiopub_clin_sample$SAMPLE_ID)
sum(is.na(posicion))



# ******************************************************************************
# genes diferencialmente expresados entre clústers HCL_P
# ******************************************************************************