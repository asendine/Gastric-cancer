# Install packages CBIO --------------------------------------------------------
bioc_packages <- c(
  "cBioPortalData",
  "TCGAbiolinks",
  "SummarizedExperiment"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# install-packages-CRAN --------------------------------------------------------
packages <- c("data.table", "dataframeexplorer", "devtools")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# libraries --------------------------------------------------------------------
library(cBioPortalData)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(pheatmap)
library(data.table)
library(dataframeexplorer)
library(devtools)
library(ComplexHeatmap)
library(cluster)
library(edgeR)
library(limma)
library(grid)
library(ggplot2)

# directory definition ---------------------------------------------------------
cbio_dir <- Sys.getenv("CBIO_DATA")
tcga_dir <- Sys.getenv("TCGA_DATA")
pub_stad_dir <- file.path(cbio_dir, "stad_tcga_pub")
gdc_stad_dir <- file.path(cbio_dir, "stad_tcga_gdc")

# Retreiving cBioPortal data 1--------------------------------------------------
cbiopub_clin_sample <- read.delim(file.path(pub_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
cbiopub_clin_pat <- read.delim(file.path(pub_stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)
cbiopub_cna <- read.delim(file.path(pub_stad_dir, "data_cna.txt"),check.names = FALSE)
cbiopub_cna_lin <- read.delim(file.path(pub_stad_dir, "data_linear_cna.txt"),
                              check.names = FALSE)
cbiopub_mrna <- read.delim(file.path(pub_stad_dir, "data_mrna_seq_v2_rsem.txt"),
                           check.names = FALSE)
cbiopub_zmrna <- read.delim(file.path(pub_stad_dir, "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"),
                            comment.char = "#", check.names = FALSE)
cbiopub_mut <- read.delim(file.path(pub_stad_dir, "data_mutations.txt"),comment.char = "#",
                          check.names = FALSE)

# Retreiving cBioPortal data 2---------------------------------------------------
cbiogdc_clin_sample <- read.delim(file.path(gdc_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
cbiogdc_clin_pat <- read.delim(file.path(gdc_stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)
cbiogdc_cna <- read.delim(file.path(gdc_stad_dir, "data_cna.txt"), check.names = FALSE)
cbiogdc_mrna_tpm <- read.delim(file.path(gdc_stad_dir, "data_mrna_seq_tpm.txt"),
                                check.names = FALSE)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# Retreiving TCGA data from .rds files -----------------------------------------
clinical_tcga <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_clinical.rds"))
cnv_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_cnv_se.rds"))
maf <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_maf.rds"))
methylation_450K_se <- readRDS(
  file.path(tcga_dir, "Prepared", "TCGA_STAD_methylation_450K_se.rds"))
mirna_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_mirna_se.rds"))
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))

# OBJETIVO 1: CLASIFICAR LOS 443 CASOS DE GDCportal
# 1 - MEDIANTE RÉPLICA ALGORITMO (DIFÍCIL?)
# 1.1 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING Platform-Specific
# 1.1.1 - (S2)SCNA (somatic copy number alterations): SCNAs vs localización cromosoma.
        # Luego cluster jerárqico mediante GISTIC2.0 https://github.com/broadinstitute/gistic2
# 1.1.2 - ... hasta (S7)
# 1.2 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING molecular data with iCluster+ (paquete bioc)

********************************************************************************
# EMPEZAMOS CON RNA-SEQ*********************************************************
********************************************************************************

# ------------------------------------------------------------------------------
# DESCUBRIMIENTO DE LOS CLÚSTERS # ---------------------------------------------
# ------------------------------------------------------------------------------

assayNames(rnaseq_se)
# unstranded: raw count
# stranded first: raw count a partir del 1a cDNA
# stranded second: raw count a partir del primer cDNA
# tpm_unstrand: Transcripts Per Million calculated using unstranded raw counts.
# fpkm_unstrand: Fragments Per Kilobase of transcript per Million mapped reads calculated from unstranded counts.
# fpkm_uq_unstrand: Upper Quartile (UQ) normalized FPKM.

tpm <- assay(rnaseq_se, "tpm_unstrand") # todas las muestras del proyecto TCGA-STAD con datos rna-seq
dim(tpm)
# esto es 448 muestras
summary(colSums(tpm)) # paso extra cortesía de GPT, debería sumar 10^6
sum(tpm == 0, na.rm = TRUE)
# hay bastantes 0
# eliminamos las muestras con valor 0 en todos sus genes
tpm_filt <- tpm[, colSums(tpm > 0) > 0, drop = FALSE]
dim(tpm_filt)
detect_dupl_cols(tpm_filt, return_type = "col_names", duplicate_col = "right")
# no se elimina nada, todo parece ok

# ahora se debería ver si todas las muestras que hay son de tumor o tejido normal
# se ve que en los nombres de las columnas, es decir muestras, la posición 14-15
# indica si es tumor (01) o tejido normal (11)
tpm_tumor <- tpm_filt[, substr(colnames(tpm_filt), 14, 15) == "01", drop = FALSE]
dim(tpm_tumor)
# hay 412 muestras de tumor

# ahora se mantendran los genes con un tpm >= 1 en el 25% de las muestras
# se calcula qué número es el 25%
n_min <- ceiling(0.25*ncol(tpm_tumor))
# ahora buscamos el dataset en el que los valores sean >= 1 en 103 columnas o más
genes_exp <- rowSums(tpm_tumor >= 1) >= n_min
# nos queda un objeto con valores TRUE/FALSE con su fila según si cumple o no las condiciones
tpm_filt <- tpm_tumor[genes_exp, , drop = FALSE]
dim(tpm_filt)
# quedan 21k genes
# aquí se deberían evaluar si hay outliers o batch effect
# se transforman los datos a log
tpm_filt_log <- log2(tpm_filt + 1)
# Ahora se filtra por variabilidad de expresión (median absolute deviation) buscando el 25% más variable
gene_mad <- apply(tpm_filt_log, 1, mad) # esto da un vector de mad aplicado a cada fila (1) de la matrix indicada
# igual que genes_exp pero numérico en vez de lógico
top <- ceiling(0.25*nrow(tpm_filt_log)) # se obtiene el nº de genes = 25%
ind <- order(gene_mad, decreasing = TRUE)[seq_len(top)] 
# order devuelve la posición
# seq_len crea un vector del 1 al 5349
# se obtienen las posiciones (filas) del 1 al 5349

# se limita el dataset a las filas que hay en ind
tpm_filt <- tpm_filt_log[ind, , drop=FALSE]
dim(tpm_filt)
# quedan 5349 genes

# se realiza un análisis del efecto batch simplemente por añadir calidad al workflow
# ya que sabemos por la documentación oficial que solamente hay un ligero batch effect en 
# datos de miRNA (1/4 analizados) pero no es trascendente. También se puede hacer la 
# visualización directa en la web de MDAnderson sin tener que hacerlo en R (tienen un
# web browser para ello, PCA Plus).
# A partir del nombre de las muestras se obtiene el id de placa y el TSS (centro)
# que estan en las posiciones 6 y 2 respectivamente (TCGA-BR-4257-01A-01R-1131-13).

barcode_parts <- do.call(rbind, strsplit(colnames(rnaseq_se), "-", fixed = TRUE))
colData(rnaseq_se)$TSS <- factor(barcode_parts[, 2])
colData(rnaseq_se)$PlateId <- factor(barcode_parts[, 6])

batch_info <- as.data.frame(colData(rnaseq_se)[,c("TSS", "PlateId")])
head(batch_info)

pca_batch <- prcomp(t(tpm_filt), center = TRUE, scale. = FALSE)
variance_explained <- 100*pca_batch$sdev^2/sum(pca_batch$sdev^2)
sample_position <- match(rownames(pca_batch$x), colnames(rnaseq_se))
pca_batch_df <- data.frame(sample = rownames(pca_batch$x),
                           PC1 = pca_batch$x[, 1],
                           PC2 = pca_batch$x[, 2],
                           PlateId = colData(rnaseq_se)$PlateId[sample_position],
                           TSS = colData(rnaseq_se)$TSS[sample_position])

ggplot(pca_batch_df, aes(x = PC1, y = PC2, color = TSS)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "TSS") +
  theme_bw()

ggplot(pca_batch_df, aes(x = PC1, y = PC2, color = PlateId)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "ID Placa") +
  theme_bw()

# 1- separar estas muestras en clústers

# como queremos clusterizar muestras (columnas) habría que transponer la matriz 
# para que la disimilitud se calcule entre los genes

# primero se transpone la matriz, se hace el z-score y se vuelve a transponer de nuevo.
tpm_z <- t(scale(t(tpm_filt)))
sum(is.na(tpm_z))

# clustering más habitual: euclidean + ward.D2. Distancia + clust
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")
plot(hcl, labels = FALSE, hang = -1, main = "Hierarchical clust: euc + wardd2", xlab = "muestras")

# clustering 2: 1-Pearson + average, mide correlaciones entre rows y hace el average?
# se calcula la correlacion y la dist
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
# se hace el cluster
hcl_p <- hclust(dist_pearson, method = "average")
plot(hcl_p, labels = FALSE, hang = -1, main = "Hierarchical clust: 1-Pearson + average", xlab = "muestras")

# aparentemente k = 4 parece lo mejor en ambas opciones
# veamos la n en cada clúster para cada método
clusters <- cutree(hcl, k = 4)
table(clusters)

clustersp <- cutree(hcl_p, k = 4)
table(clustersp)

# se forman 4 clusters escalonados en ambos métodos
# por ahora se usará eucliden+ward.D2
# ------------------------------------------------------------------------------
# CARACTERIZACIÓN DE LOS CLÚSTERS # --------------------------------------------
# ------------------------------------------------------------------------------
# 2- buscar los genes diferencialmente más expresados de cada clúster (~10)
# ahora hay que hacer el análisis de expresión diferencial por cluster
# obtenemos los datos de conteos crudos para poder usar mejor edger y limma
counts <- assay(rnaseq_se, "unstranded")
dim(counts)
# se obtienen las muestras de tumor primario
counts <- counts[, names(clusters), drop = FALSE]
dim(counts)
stopifnot(identical(colnames(counts), names(clusters))) # para ver facilmente que los nombres sean iguales
# se crea el objeto DGEList
# primero se obtiene el vector categórico que hace de "meta" de las muestras
cluster <- factor(clusters[colnames(counts)], levels = 1:4, labels = c("C1","C2","C3","C4"))
# esta línea es compleja... te dice que la muestra x pertenece al clúster y y así con todas
# se crea el DGEList
dge <- DGEList(counts = counts, group = cluster)
# se hace un filtrado de genes preventivo con filterByExpr(). Filtrado por defecto, sin más.
keep <- filterByExpr(dge, group = cluster)
# si se elimina un gen, recalcula los counts por muestra, esto es cosa de GPT, nunca lo he usado:
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
# se hace la matriz de diseño siguiendo el modelo cluster
design <- model.matrix(~ 0 + cluster)
head(design)
# se cambian los nombres de las columnas
colnames(design) <- levels(cluster)
head(design)
# se crea el voom para hacer el modelo lineal
v <- voom(dge, design)
# se ajusta el modelo lineal
fit <- lmFit(v, design)
# se organizan los contrastes entre clusters, se compara la expresión de cada clúster
# con la media de expresiones del resto... es una forma pero hay otras, tipo ver
# genes expresados diferencialmente exclusivos de cada cluster, podria ser interesante?
# sintaxis es diferentes contrastes como un vector + matriz de diseño con los parámetros como columnas
contrasts <- makeContrasts(C1all = C1 - (C2+C3+C4)/3,
                           C2all = C2 - (C1+C3+C4)/3,
                           C3all = C3 - (C1+C2+C4)/3,
                           C4all = C4 - (C1+C2+C3)/3, levels = design)
# contrasts.fit() coge un ajuste de un modelo lineal (fit) y una matriz con filas
# que equivalen a las columnas de los coeficientes de fit y columnas que equivalen a los contrastes?
# realiza el ajuste con contrastes
fit_clusters <- contrasts.fit(fit, contrasts)
# treat testea cambios diferentes a un límite dado (lfc > x), eBayes() testea si hay diferencias respecto a 0
fit_treat <- treat(fit_clusters, lfc = 0.95)
# se obtienen los top-ranked genes más expresados diferencialmente de cada contraste
# lfc = 1 ya filtra por valor absoluto, "logFC" también
genes_C1 <- topTreat(fit_treat, coef = "C1all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C2 <- topTreat(fit_treat, coef = "C2all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C3 <- topTreat(fit_treat, coef = "C3all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C4 <- topTreat(fit_treat, coef = "C4all", number = Inf, p.value = 0.05, sort.by = "logFC")
info_C1 <- head(genes_C1, 10) # si lfc = 0 salen muchos genes, pero si lfc = 1 solo se filtran 5!!
info_C2 <- head(genes_C2, 10) # finalmente indico 0.96 para que salgan 11 genes en C1
info_C3 <- head(genes_C3, 10) # con esto tenemos los genes en formato gene_id
info_C4 <- head(genes_C4, 10)

# 3- hacer el heatmap
# Primero se prepara la información principal de cada clúster
gene_info <- as.data.frame(rowData(rnaseq_se))
head(gene_info)
# C1 ----------------------------------------------------------
info_C1$gene_id <- rownames(info_C1)
# la siguiente línea con match busca el valor de gene_name en gene_info antes definido
# que corresponde entre el gene_id de info_C1 (ENSG...) y las filas de gene_info (ENSG...)
info_C1$gene_symbol <- gene_info$gene_name[match(info_C1$gene_id, rownames(gene_info))]
info_C1$gene_cluster <- "C1"
head(info_C1)
# C2 ----------------------------------------------------------
info_C2$gene_id <- rownames(info_C2)
info_C2$gene_symbol <- gene_info$gene_name[match(info_C2$gene_id, rownames(gene_info))]
info_C2$gene_cluster <- "C2"
head(info_C2)
# C3 ----------------------------------------------------------
info_C3$gene_id <- rownames(info_C3)
info_C3$gene_symbol <- gene_info$gene_name[match(info_C3$gene_id, rownames(gene_info))]
info_C3$gene_cluster <- "C3"
head(info_C3)
# C4 ----------------------------------------------------------
info_C4$gene_id <- rownames(info_C4)
info_C4$gene_symbol <- gene_info$gene_name[match(info_C4$gene_id, rownames(gene_info))]
info_C4$gene_cluster <- "C4"
head(info_C4)
# nos quedamos solo con las columnas que queremos y lo unimos todo
columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")
info_C1 <- info_C1[, columns_keep]
info_C2 <- info_C2[, columns_keep]
info_C3 <- info_C3[, columns_keep]
info_C4 <- info_C4[, columns_keep]
genes_heatmap_info <- rbind(info_C1, info_C2, info_C3, info_C4)
dim(genes_heatmap_info)

# Segundo se construye la matriz de expresión
# resulta que el objeto "v" creado por voom contiene la relación entre el gene_id
# y las muestras. Los valores no sé exactamente qué son
gene_ids <- genes_heatmap_info$gene_id
# se filtran los 40 genes que nos interesan
mat_heatmap <- v$E[gene_ids, , drop = FALSE]
# y se calcula el z-score de cada valor
mat_z <- mat_heatmap

# básicamente para cada elemento en mat_heatmap se calcula la media y la sd
# entonces se calcula el z-score de cada valor con la media y la sd y se sustituye
# en mat_z:
for (i in 1:nrow(mat_heatmap)) {
  gene_mean <- mean(mat_heatmap[i, ])
  gene_sd <- sd(mat_heatmap[i, ])
  mat_z[i, ] <- (mat_heatmap[i, ] - gene_mean)/gene_sd
}

# por último se crea la agrupación por clusters para las filas y las columnas
gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = c("C1", "C2", "C3", "C4"))

sample_cluster <- factor(cluster[colnames(mat_z)], levels = c("C1", "C2", "C3", "C4"))


# heatmap
# modificado para obtener los clusters ordenados y añadidos debajo como bloques
# por colores
cluster_colors <- c(
  "C1" = "#E64B35",
  "C2" = "#4DBBD5",
  "C3" = "#00A087",
  "C4" = "#F39B7F")

annotation_clusters <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors),
    labels = c("C1", "C2", "C3", "C4"),
    labels_gp = gpar(col = "white", fontface = "bold")
  ),
  which = "column"
)
# h1
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  column_split = sample_cluster,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  bottom_annotation = annotation_clusters,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm"),
  row_title = "Genes DE")

# finalmente se relaciona con los molecular subtypes obtenidos en las 290 y pico muestras

sample_id_mat <- substr(colnames(mat_z), 1, 15)
posicion <- match(sample_id_mat, cbiopub_clin_sample$SAMPLE_ID)
sum(is.na(posicion))
# hay 138 NAs, es decir 274 muestras etiquetadas con subtype de las 412, pero en 
# cbio constan 295 muestras etiquetadas. Hay 21 muestras que no estan entre las 412?

molecular_subtype <- cbiopub_clin_sample$MOLECULAR_SUBTYPE[posicion]
names(molecular_subtype) <- colnames(mat_z)
head(molecular_subtype)
subtype_colors <- c(
  "CIN" = "#3C5488",
  "EBV" = "#00A087",
  "GS" = "#FFD700",
  "MSI" = "#E64B35"
)

# h2
annotation_bottom <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[c("C1", "C2", "C3", "C4")]),
    labels = c("C1", "C2", "C3", "C4"),
    labels_gp = gpar(col = "white",fontface = "bold")),
  Molecular_subtype = molecular_subtype,
  col = list(Molecular_subtype = subtype_colors),
  na_col = "#D9D9D9",
  annotation_height = unit(c(6, 4),"mm"),
  show_annotation_name = c(Cluster = FALSE, 
                           Molecular_subtype = TRUE),
  which = "column")

Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  column_split = sample_cluster,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  bottom_annotation = annotation_bottom,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)
