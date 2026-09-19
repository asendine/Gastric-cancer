# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'readxl', 'here', 'pROC',
               'tidyestimate', 'tidyr', 'uwot')
for (i in libraries) {
  library(i, character.only = TRUE)
}

# ******************************************************************************
# DIRECTORY SETUP
# ******************************************************************************
tcga_dir <- Sys.getenv("TCGA_DATA")
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))
tpm_ind <- readRDS(here("misc", "data", "tpm_ind.rds"))
# ******************************************************************************

# ******************************************************************************
# BACTH EFFECT
# ******************************************************************************
# A partir del nombre de las muestras se obtiene el id de placa y el TSS (origen muestra)
# que estan en las posiciones 6 y 2 respectivamente (TCGA-BR-4257-01A-01R-1131-13).
barcode_parts <- do.call(rbind, strsplit(colnames(rnaseq_se), "-", fixed = TRUE))
colData(rnaseq_se)$TSS <- factor(barcode_parts[, 2])
colData(rnaseq_se)$PlateId <- factor(barcode_parts[, 6])
batch_info <- as.data.frame(colData(rnaseq_se)[,c("TSS", "PlateId")])

pca_batch <- prcomp(t(tpm_ind), center = TRUE, scale. = FALSE)
variance_explained <- 100*pca_batch$sdev^2/sum(pca_batch$sdev^2)
sample_position <- match(rownames(pca_batch$x), colnames(rnaseq_se))
pca_batch_df <- data.frame(sample = rownames(pca_batch$x),
                           PC1 = pca_batch$x[, 1],
                           PC2 = pca_batch$x[, 2],
                           PlateId = batch_info$PlateId[sample_position],
                           TSS = batch_info$TSS[sample_position])

pca_background <- pca_batch_df[, c("PC1", "PC2")]

ggplot(pca_batch_df, aes(PC1, PC2)) + 
  geom_point(data = pca_background, aes(PC1, PC2), inherit.aes = FALSE, 
             color = "grey85", size = 0.7) + 
  geom_point(color = "#D55E00", size = 1.3, alpha = 0.9) + 
  facet_wrap(~ TSS, ncol = 5) +
  labs(title = "PCA: distribución de las muestras por TSS", 
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggplot(pca_batch_df, aes(PC1, PC2)) + 
  geom_point(data = pca_background, aes(PC1, PC2), inherit.aes = FALSE, 
             color = "grey85", size = 0.7) + 
  geom_point(color = "#0072B2", size = 1.3, alpha = 0.9) + 
  facet_wrap(~ PlateId, ncol = 5) +
  labs(title = "PCA: distribución de las muestras por ID de placa", 
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# ******************************************************************************
# OUTLIERS
# ******************************************************************************
pca_outliers <- pca(mat = tpm_ind, center = TRUE, 
                              scale = FALSE, removeVar = NULL)
screeplot(pcaobj = pca_outliers, components = getComponents(pca_outliers, 1:10), 
          title = "PCA mRNA")
pairsplot(pcaobj = pca_outliers, components = getComponents(pca_outliers, 1:3),
          triangle  = TRUE)

# ******************************************************************************
# GENES DOMINANTES
# ******************************************************************************
gene_info <- as.data.frame(rowData(rnaseq_se))
gene_index <- match(rownames(tpm_ind), rownames(rnaseq_se))
gene_map <- data.frame(ensembl_id = rownames(tpm_ind),
                       gene_name = gene_info$gene_name[gene_index])

# dado un pca, una componente, la leyenda y el nº de genes a comprobar
get_top_genes <- function(pca_object, pc, gene_map, n = 10) { 
  loading_values <- pca_object$loadings[, pc] 
  top_idx <- order(abs(loading_values), decreasing = TRUE)[seq_len(n)]
  ensembl_ids <- rownames(pca_object$loadings)[top_idx]
  data.frame(PC = pc,
             gene_name = gene_map$gene_name[match(ensembl_ids, gene_map$ensembl_id)],
             ensembl_id = ensembl_ids,
             loading = loading_values[top_idx],
             contribution_pct = 100*loading_values[top_idx]^2)
}

cont_PC1 <- get_top_genes(pca_outliers, "PC1", gene_map, n = 10)
cont_PC2 <- get_top_genes(pca_outliers, "PC2", gene_map, n = 10)
cont_PC1
cont_PC2

# umap
idx <- match(substr(colnames(tpm_ind), 1, 15), cbiopub_clin_sample$SAMPLE_ID)
set.seed(1234)
umap_result <- umap(X = t(as.matrix(tpm_ind)), 
                    n_neighbors = 30, 
                    min_dist    = 0.1, # cuánto pueden compactarse los puntos en el mapa.
                    metric      = "correlation", 
                    verbose     = FALSE)
umap_df <- data.frame(UMAP1   = umap_result[, 1],
                      UMAP2   = umap_result[, 2],
                      subtipo = cbiopub_clin_sample$MOLECULAR_SUBTYPE[idx])

umap_df$subtipo <- as.character(umap_df$subtipo)
umap_df$subtipo[is.na(umap_df$subtipo)] <- "Sin etiqueta"

ggplot(umap_df, aes(UMAP1, UMAP2, colour = subtipo)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = c(
    CIN = "darkorchid4",
    EBV = "#00A087",
    GS  = "#B8860B",
    MSI = "#FFFF00",
    "Sin etiqueta" = "grey75")) +
  coord_equal() +
  labs(title = "UMAP TCGA-STAD", x = "UMAP 1", y = "UMAP 2", 
       colour = "Subtipo TCGA") +
  theme_classic()
