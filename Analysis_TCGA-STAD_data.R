# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'writexl', 'uwot')
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
tpm_ind <- readRDS(here("misc", "data", "tpm_ind.rds"))
tpm_z <- t(scale(t(tpm_ind)))
nmf_rank <- readRDS("misc/data/nmf_rank.rds")
cc_input <- as.matrix(tpm_z)
cc_results <- ConsensusClusterPlus(
  d             = cc_input,
  maxK          = 6,            # evalúa k = 2,...,6
  reps          = 500,           # 500 remuestreos
  pItem         = 0.80,          # 80% de las muestras en cada repetición
  pFeature      = 0.80,          # 80% de los genes en cada repetición
  clusterAlg    = "hc",          
  distance      = "euclidean",
  innerLinkage  = "ward.D2",     # algoritmo aplicado en cada repetición
  finalLinkage  = "average",     # agrupación final de la matriz de consenso
  seed          = 1234,
  plot          = NULL,
  writeTable    = FALSE,
  verbose       = TRUE)

cc_resultsp <- ConsensusClusterPlus(
  d             = cc_input,
  maxK          = 6,
  reps          = 500,           
  pItem         = 0.80,          
  pFeature      = 0.80,          
  clusterAlg    = "hc",          
  distance      = "pearson",
  innerLinkage  = "average",     
  finalLinkage  = "average",     
  seed          = 1234,
  plot          = NULL,
  writeTable    = FALSE,
  verbose       = TRUE)

# ******************************************************************************
# CARACTERIZACIÓN DE LOS CLÚSTERS
# ******************************************************************************
# el siguiente paso es encontrar las "core samples", aquellas con anchura de silueta
# más positiva ~mayor similaridad con su cluster.
# "Samples most representative of the clusters, hereby called core samples were 
# identified based on positive silhouette width, indicating higher similarity to 
# their own class than to any other class member. Core samples were used to select 
# differentially expressed marker genes for each subtype by comparing the subclass 
# versus the other subclasses, using Student's t-test." La realidad es que usando 
# cualquier valor de k se han obtenido valores extremadamente buenos de anchura de
# silueta. Aquí surge otro dilema, plantear un criterio arbitrario para filtrar muestras
# o dejarlas todas y ya está... no me parece mal en un contexto exploratorio
# obtener solamente las muestras que más definan un cluster sin pasarnos. El criterio
# silueta > 0 no discrimina nada, empezaría a ser interesante con sil > 0.6, en el que
# se eliminarían algunas muestras, quizá un 30%, que ya está bien.
k <- 4
fit_nmf <- nmf_rank$fit[[as.character(k)]]
asignaciones <- list(NMF = predict(fit_nmf, what = "consensus"))

# Para incluir ConsensusClusterPlus, descomentar estas dos lineas:
asignaciones$CC_euclidean <- cc_results[[k]]$consensusClass
asignaciones$CC_pearson <- cc_resultsp[[k]]$consensusClass

titulos <- c(NMF = "NMF: consenso",
             CC_euclidean = "ConsensusClusterPlus: euclidean",
             CC_pearson = "ConsensusClusterPlus: Pearson")
con <- paste0("C", seq_len(k))
counts_all <- assay(rnaseq_se, "unstranded")
gene_info <- as.data.frame(rowData(rnaseq_se))

# Cada cluster frente a la media (no ponderada) de los restantes.
mc <- vapply(seq_len(k), function(i) {
  paste0(con[i], " - (", paste(con[-i], collapse = " + "), ") / ", k - 1)
}, character(1))

cluster_colors <- setNames(hcl.colors(k, palette = "Dark 3"), con)
subtype_colors <- c(CIN = "darkorchid4", EBV = "#00A087",
                    GS = "#B8860B", MSI = "#FFFF00")

resultados_extra <- list()
heatmaps_extra <- list()

for (metodo in names(asignaciones)) {
  grupos <- asignaciones[[metodo]]
  stopifnot(!is.null(names(grupos)), all(names(grupos) %in% colnames(counts_all)),
            !anyNA(grupos), setequal(as.character(grupos), as.character(seq_len(k))))
  counts <- counts_all[, names(grupos), drop = FALSE]
  cluster <- factor(grupos, levels = seq_len(k), labels = con)
  names(cluster) <- names(grupos)
  stopifnot(identical(colnames(counts), names(cluster)))
  
  dge <- DGEList(counts = counts, group = cluster)
  keep <- filterByExpr(dge, group = cluster)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ 0 + cluster)
  colnames(design) <- con
  v <- voom(dge, design)
  fit_lm <- lmFit(v, design)
  
  contrasts_mat <- makeContrasts(contrasts = mc, levels = design)
  colnames(contrasts_mat) <- paste0(con, "all")
  fit_ebayes <- eBayes(contrasts.fit(fit_lm, contrasts_mat))
  
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
  
  annotation_ms <- HeatmapAnnotation(
    Cluster = anno_block(gp = gpar(fill = cluster_colors), labels = con, 
                         labels_gp = gpar(col = "white", fontface = "bold")),
    Molecular_subtype = molecular_subtype,
    col = list(Molecular_subtype = subtype_colors),
    na_col = "white",
    annotation_height = unit(c(6, 4), "mm"),
    show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE))
  
  heatmaps_extra[[metodo]] <- Heatmap(
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
  
  resultados_extra[[metodo]] <- list(cluster = cluster, genes_Ck = genes_Ck,
                                     genes_heatmap_info = genes_heatmap_info,
                                     fit_ebayes = fit_ebayes)
}

# Heatmap NMF; las otras llamadas se activan si se incluyeron sus asignaciones.
draw(heatmaps_extra$NMF)
draw(heatmaps_extra$CC_euclidean)
draw(heatmaps_extra$CC_pearson)

