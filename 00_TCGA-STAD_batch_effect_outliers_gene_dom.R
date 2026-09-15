# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'readxl', 'here', 'pROC',
               'tidyestimate', 'tidyr')
for (i in libraries) {
  library(i, character.only = TRUE)
}

# ******************************************************************************
# DIRECTORY SETUP
# ******************************************************************************
tcga_dir <- Sys.getenv("TCGA_DATA")
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))
# ******************************************************************************


# ******************************************************************************
# BACTH EFFECT
# ******************************************************************************
# Publicación TCGA describe un ligero batch effect en datos de miRNA (1/4 analizados)
# pero no es trascendente. También se puede hacer la visualización directa en la 
# web de MDAnderson sin tener que hacerlo en R (tienen un web browser para ello, PCA Plus).
# ******************************************************************************
# A partir del nombre de las muestras se obtiene el id de placa y el TSS (origen muestra)
# que estan en las posiciones 6 y 2 respectivamente (TCGA-BR-4257-01A-01R-1131-13).
barcode_parts <- do.call(rbind, strsplit(colnames(rnaseq_se), "-", fixed = TRUE))
colData(rnaseq_se)$TSS <- factor(barcode_parts[, 2])
colData(rnaseq_se)$PlateId <- factor(barcode_parts[, 6])
batch_info <- as.data.frame(colData(rnaseq_se)[,c("TSS", "PlateId")])

pca_batch <- prcomp(t(tpm_filt_log), center = TRUE, scale. = TRUE)
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
# para buscar outliers se realiza un PCA. Se usa PCAtools
pca_outliers <- PCAtools::pca(mat = tpm_filt_log, center = TRUE, 
                              scale = FALSE, removeVar = NULL)

# screeplot para ver % var explicada por cada componente
PCAtools::screeplot(pcaobj = pca_outliers, 
                    components = PCAtools::getComponents(pca_outliers, 1:20), 
                    title = "PCA mRNA")

# no está mal hacer un pairsplot entre las diferentes componentes, pero si hay
# muchas muestras y queremos ver muchas componentes, entonces no se verá bien
PCAtools::pairsplot(pcaobj = pca_outliers, 
                    components = PCAtools::getComponents(pca_outliers, 1:4),
                    triangle  = TRUE) # a partir de 4 la cosa empeora
# a partir de los plots anteriores, no se aprecian outliers


# ******************************************************************************
# GENES DOMINANTES
# ******************************************************************************
# ahora se puede buscar si hay genes que dominen alguna componente concreta
# se obtiene la info de genes
gene_info <- as.data.frame(rowData(rnaseq_se))
# se obtiene la posición de las filas de los genes que nos interesan con match
gene_index <- match(rownames(tpm_filt_log), rownames(rnaseq_se))
# se crea el dataset específico con ambas notaciones
gene_map <- data.frame(ensembl_id = rownames(tpm_filt_log),
                       gene_name = gene_info$gene_name[gene_index])
# a partir de aquí se crea una función para determinar los genes que más contribuyen
# dado un pca, una componente, la leyenda y el nº de genes a comprobar
get_top_genes <- function(pca_object, pc, gene_map, n = 10) { 
  loading_values <- pca_object$loadings[, pc] 
  # se obtienen los pesos
  top_idx <- order(abs(loading_values), decreasing = TRUE)[seq_len(n)]
  # se obtienen los n pesos por valor absoluto decreciente
  ensembl_ids <- rownames(pca_object$loadings)[top_idx]
  # se obtienen los ids de los genes con esos pesos
  data.frame(PC = pc,
             gene_name = gene_map$gene_name[match(ensembl_ids, gene_map$ensembl_id)],
             ensembl_id = ensembl_ids,
             loading = loading_values[top_idx],
             contribution_pct = 100*loading_values[top_idx]^2)
  # se hace un dataframe con la pc elegida, se busca el gene name con la leyenda indicada,
  # se indica el ensemble id, el peso obtenido y las contribuciones
}

# aquí se puede ver qué genes dominan las componentes en este caso la 1 y 2
cont_PC1 <- get_top_genes(pca_outliers, "PC1", gene_map, n = 10)
cont_PC2 <- get_top_genes(pca_outliers, "PC2", gene_map, n = 10)
cont_PC1
cont_PC2

# pendiente ver si estos genes se asocian con alguna característica concreta