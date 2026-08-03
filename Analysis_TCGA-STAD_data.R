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

# por último se crea un "groups" con los clusters
gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = c("C1", "C2", "C3", "C4"))

sample_cluster <- cluster[colnames(mat_z)]

# el heatmap
Heatmap(
  mat_z,
  name = "Z-score",
  
  row_labels = genes_heatmap_info$gene_symbol,
  row_split = gene_cluster,
  column_split = sample_cluster,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  show_column_names = FALSE,
  row_title = "Genes DE"
)
