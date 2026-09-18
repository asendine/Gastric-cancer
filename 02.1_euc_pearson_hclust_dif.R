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


# Retreiving cBioPortal data****************************************************
cbiopub_clin_sample <- read.delim(file.path(pub_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
cbiopub_clin_pat <- read.delim(file.path(pub_stad_dir, "data_clinical_patient.txt"),
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

k <- 4

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

coefs <- character(k)
genes_Ck <- vector("list", k)
info_Ck <- vector("list", k)
names(genes_Ck) <- paste0("C", seq_len(k))
names(info_Ck) <- names(genes_Ck)

gene_info <- as.data.frame(rowData(rnaseq_se))
columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")

for (n in seq_len(k)) {
  coefs[n] <- paste0("C", n, "all")
  
  genes_Ck[[n]] <- topTable(fit_ebayes, coef = coefs[n], number = Inf, 
                            p.value = 0.05, sort.by = "logFC")
  
  info_Ck[[n]] <- head(genes_Ck[[n]], 10)
  info_Ck[[n]]$gene_id <- rownames(info_Ck[[n]])
  info_Ck[[n]]$gene_symbol <- gene_info$gene_name[match(info_Ck[[n]]$gene_id, 
                                                        rownames(gene_info))]
  info_Ck[[n]]$gene_cluster <- rep(paste0("C", n), nrow(info_Ck[[n]]))
  info_Ck[[n]] <- info_Ck[[n]][, columns_keep]
}

# pasa cada tabla de la lista info_Ck a rbind()
genes_heatmap_info <- do.call(rbind, info_Ck)
gene_ids <- genes_heatmap_info$gene_id

mat_heatmap <- v$E[gene_ids, , drop = FALSE]


# ******************************************************************************
# heatmap
# ******************************************************************************

mat_z <- mat_heatmap
for (i in 1:nrow(mat_heatmap)) {
  gene_mean <- mean(mat_heatmap[i, ])
  gene_sd <- sd(mat_heatmap[i, ])
  mat_z[i, ] <- (mat_heatmap[i, ] - gene_mean)/gene_sd
}

anyDuplicated(genes_heatmap_info$gene_id)

gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = paste0("C", seq_len(k)))
cluster_colors <- setNames(hcl.colors(k, palette = "Dark 3"), paste0("C", seq_len(k)))
row_annotation <- rowAnnotation(Clusters = gene_cluster,
                                col = list(Clusters = cluster_colors))
sample_cluster <- factor(cluster[colnames(mat_z)],
                         levels = paste0("C", seq_len(k)))

sample_id_mat <- substr(colnames(mat_z), 1, 15)
posicion <- match(sample_id_mat, cbiopub_clin_sample$SAMPLE_ID)
sum(is.na(posicion))

molecular_subtype <- cbiopub_clin_sample$MOLECULAR_SUBTYPE[posicion]
names(molecular_subtype) <- colnames(mat_z)
head(molecular_subtype)
subtype_colors <- c("CIN" = "darkorchid4", "EBV" = "#00A087", "GS" = "#B8860B", "MSI" = "#FF69B4")

# anotación inferior 2: Molecular subtype
annotation_ms <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[c(paste0("C", 1:k))]),
    labels = c(paste0("C", 1:k)),
    labels_gp = gpar(col = "white",fontface = "bold")),
  Molecular_subtype = molecular_subtype,
  col = list(Molecular_subtype = subtype_colors),
  na_col = "white",
  annotation_height = unit(c(6, 4),"mm"),
  show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE),
  which = "column")

# cluster genes (separacion genes por tipo de cluster)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  row_split = gene_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)

# row clustering 1-pearson + average (como se expresan en general)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  row_dend_reorder = FALSE,
  row_split = gene_cluster,
  left_annotation = row_annotation,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm"),
  right_annotation = row_annotation,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "average"
)


# ******************************************************************************
# genes diferencialmente expresados entre clústers HCL_P
# ******************************************************************************

counts <- assay(rnaseq_se, "unstranded")
counts <- counts[, names(cluster_hcl_p), drop = FALSE]
stopifnot(identical(colnames(counts), names(cluster_hcl_p)))

cluster <- factor(cluster_hcl_p[colnames(counts)], levels = 1:k, labels = paste0("C", 1:k))
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

coefs <- character(k)
genes_Ck <- vector("list", k)
info_Ck <- vector("list", k)
names(genes_Ck) <- paste0("C", seq_len(k))
names(info_Ck) <- names(genes_Ck)

gene_info <- as.data.frame(rowData(rnaseq_se))
columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")

for (n in seq_len(k)) {
  coefs[n] <- paste0("C", n, "all")
  
  genes_Ck[[n]] <- topTable(fit_ebayes, coef = coefs[n], number = Inf, 
                            p.value = 0.05, sort.by = "logFC")
  
  info_Ck[[n]] <- head(genes_Ck[[n]], 10)
  info_Ck[[n]]$gene_id <- rownames(info_Ck[[n]])
  info_Ck[[n]]$gene_symbol <- gene_info$gene_name[match(info_Ck[[n]]$gene_id, 
                                                        rownames(gene_info))]
  info_Ck[[n]]$gene_cluster <- rep(paste0("C", n), nrow(info_Ck[[n]]))
  info_Ck[[n]] <- info_Ck[[n]][, columns_keep]
}

# pasa cada tabla de la lista info_Ck a rbind()
genes_heatmap_info <- do.call(rbind, info_Ck)
gene_ids <- genes_heatmap_info$gene_id

mat_heatmap <- v$E[gene_ids, , drop = FALSE]


# ******************************************************************************
# heatmap
# ******************************************************************************

mat_z <- mat_heatmap
for (i in 1:nrow(mat_heatmap)) {
  gene_mean <- mean(mat_heatmap[i, ])
  gene_sd <- sd(mat_heatmap[i, ])
  mat_z[i, ] <- (mat_heatmap[i, ] - gene_mean)/gene_sd
}

anyDuplicated(genes_heatmap_info$gene_id)

gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = paste0("C", seq_len(k)))
cluster_colors <- setNames(hcl.colors(k, palette = "Dark 3"), paste0("C", seq_len(k)))
row_annotation <- rowAnnotation(Clusters = gene_cluster,
                                col = list(Clusters = cluster_colors))
sample_cluster <- factor(cluster[colnames(mat_z)],
                         levels = paste0("C", seq_len(k)))

sample_id_mat <- substr(colnames(mat_z), 1, 15)
posicion <- match(sample_id_mat, cbiopub_clin_sample$SAMPLE_ID)
sum(is.na(posicion))

molecular_subtype <- cbiopub_clin_sample$MOLECULAR_SUBTYPE[posicion]
names(molecular_subtype) <- colnames(mat_z)
head(molecular_subtype)
subtype_colors <- c("CIN" = "darkorchid4", "EBV" = "#00A087", "GS" = "#B8860B", "MSI" = "#FF69B4")

# anotación inferior 2: Molecular subtype
annotation_ms <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[c(paste0("C", 1:k))]),
    labels = c(paste0("C", 1:k)),
    labels_gp = gpar(col = "white",fontface = "bold")),
  Molecular_subtype = molecular_subtype,
  col = list(Molecular_subtype = subtype_colors),
  na_col = "white",
  annotation_height = unit(c(6, 4),"mm"),
  show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE),
  which = "column")

# cluster genes (separacion genes por tipo de cluster)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  row_split = gene_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)

# row clustering 1-pearson + average (como se expresan en general)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  row_dend_reorder = FALSE,
  row_split = gene_cluster,
  left_annotation = row_annotation,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm"),
  right_annotation = row_annotation,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "average"
)

# ******************************************************************************
# Diferencias hcl vs hcl_p
# ******************************************************************************
# revisar los heatmaps ya que salen iguales, revisar como se clusterizan las columnas