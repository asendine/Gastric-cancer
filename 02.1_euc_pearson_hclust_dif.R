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
tpm_ind <- readRDS(here("misc", "data", "tpm_ind.rds"))

tpm_z <- t(scale(t(tpm_ind)))
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
hcl_p <- hclust(dist_pearson, method = "average")

# ******************************************************************************
# corte - asociación de muestras a clústers
# ******************************************************************************
k <- 4

metodos <- list(euclidean = hcl, pearson = hcl_p)
titulos <- c(euclidean = "Euclidea + Ward.D2", pearson = "1-Pearson + average")
con <- paste0("C", seq_len(k))
counts_all <- assay(rnaseq_se, "unstranded")
gene_info <- as.data.frame(rowData(rnaseq_se))

# Cada cluster frente a la media (no ponderada) de los restantes.
mc <- vapply(seq_len(k), function(i) {
  paste0(con[i], " - (", paste(con[-i], collapse = " + "), ") / ", k - 1)
}, character(1))

cluster_colors <- setNames(hcl.colors(k, palette = "Dark 3"), con)
subtype_colors <- c(CIN = "darkorchid4", EBV = "#00A087", GS = "#B8860B", 
                    MSI = "#FFFF00")

resultados <- list()
heatmaps <- list()

for (metodo in names(metodos)) {
  grupos <- cutree(metodos[[metodo]], k = k)
  stopifnot(!is.null(names(grupos)), all(names(grupos) %in% colnames(counts_all)))
  counts <- counts_all[, names(grupos), drop = FALSE]
  cluster <- factor(grupos, levels = seq_len(k), labels = con)
  names(cluster) <- names(grupos)
  stopifnot(identical(colnames(counts), names(cluster)))
  
  # Expresion diferencial: filtrado, TMM, voom y limma.
  dge <- DGEList(counts = counts, group = cluster)
  keep <- filterByExpr(dge, group = cluster)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ 0 + cluster)
  colnames(design) <- con
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  contrasts_mat <- makeContrasts(contrasts = mc, levels = design)
  colnames(contrasts_mat) <- paste0(con, "all")
  fit_ebayes <- eBayes(contrasts.fit(fit, contrasts_mat))
  
  genes_Ck <- setNames(vector("list", k), con)
  info_Ck <- setNames(vector("list", k), con)
  
  for (n in seq_len(k)) {
    genes_Ck[[n]] <- topTable(fit_ebayes, coef = paste0(con[n], "all"), 
                              number = Inf, p.value = 0.05, adjust.method = "BH", 
                              sort.by = "logFC")
    info <- head(genes_Ck[[n]], 10)
    info$gene_id <- rownames(info)
    info$gene_symbol <- gene_info$gene_name[match(info$gene_id, rownames(gene_info))]
    info$gene_cluster <- rep(con[n], nrow(info))
    info_Ck[[n]] <- info[, c("gene_id", "gene_symbol", "gene_cluster", "logFC", 
                             "adj.P.Val")]
  }
  
  genes_heatmap_info <- do.call(rbind, info_Ck)
  mat_heatmap <- v$E[genes_heatmap_info$gene_id, , drop = FALSE]
  mat_z <- t(scale(t(mat_heatmap)))
  gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = con)
  sample_cluster <- cluster[colnames(mat_z)]
  posicion <- match(substr(colnames(mat_z), 1, 15), cbiopub_clin_sample$SAMPLE_ID)
  molecular_subtype <- cbiopub_clin_sample$MOLECULAR_SUBTYPE[posicion]
  annotation_ms <- HeatmapAnnotation(Cluster = anno_block(
    gp = gpar(fill = cluster_colors), labels = con, labels_gp = gpar(
      col = "white", fontface = "bold")),
    Molecular_subtype = molecular_subtype,
    col = list(Molecular_subtype = subtype_colors),
    na_col = "white",
    annotation_height = unit(c(6, 4), "mm"),
    show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE))
  
  heatmaps[[metodo]] <- Heatmap(
    mat_z,
    name = "Z-score",
    column_title = paste0(titulos[[metodo]]),
    row_labels = genes_heatmap_info$gene_symbol,
    row_names_gp = gpar(fontsize = 7),
    column_split = sample_cluster,
    row_split = gene_cluster,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    row_dend_reorder = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "average",
    left_annotation = rowAnnotation(Clusters = gene_cluster, 
                                    col = list(Clusters = cluster_colors)),
    bottom_annotation = annotation_ms,
    show_column_names = FALSE,
    column_gap = unit(2, "mm"),
    row_gap = unit(2, "mm"))
  
  resultados[[metodo]] <- list(cluster = cluster, genes_Ck = genes_Ck,
                               genes_heatmap_info = genes_heatmap_info,
                               fit_ebayes = fit_ebayes)
}

draw(heatmaps$euclidean)
draw(heatmaps$pearson)

# ******************************************************************************
# Diferencias hcl vs hcl_p
# ******************************************************************************
