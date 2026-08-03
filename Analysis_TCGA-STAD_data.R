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
  file.path(tcga_dir, "Prepared", "TCGA_STAD_methylation_450K_se.rds")
)

mirna_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_mirna_se.rds"))

rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))

# Retreiving TCGA GDC data to classify it into the different molecular subtypes


# OBJETIVO 1: CLASIFICAR LOS 443 CASOS DE GDC FROM LOS 295 DE PUB
# 1 - MEDIANTE RÉPLICA ALGORITMO (DIFÍCIL?)
# 1.1 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING Platform-Specific
# 1.1.1 - (S2)SCNA (somatic copy number alterations): SCNAs vs localización cromosoma.
        # Luego cluster jerárqico mediante GISTIC2.0 https://github.com/broadinstitute/gistic2
# 1.1.2 - ... hasta (S7)
# 1.2 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING molecular data with iCluster+ (paquete bioc)

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

# las filas son los elementos finales del clúster, como queremos clusterizar muestras
# (columnas) habría que transponer la matriz para que la disimilitud se calcule entre los genes

# primero se transpone la matriz, se hace el z-score y se vuelve a transponer de nuevo.
tpm_z <- t(scale(t(tpm_filt_log)))
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

# aparentemente k = 4 parece lo mejor
# veamos la n en cada clúster para cada método
clusters <- cutree(hcl, k = 4)
table(clusters)

clustersp <- cutree(hcl_p, k = 4)
table(clustersp)

# aquí se ve que el clúster 1 varía un poco entre los diferentes métodos. 
# me surge una duda, cómo se puede estudiar esta diferencia??? -----------------
# por ahora se usará eucliden+ward.D2

# 2- buscar los genes diferencialmente más expresados de cada clúster (~10)
# ahora hay que hacer el análisis de expresión diferencial por cluster
# obtenemos los datos de conteos crudos para poder usar mejor edger y limma
counts <- assay(rnaseq_se, "unstranded")
# se obtienen solo las muestras de clusters
counts <- counts[, names(clusters), drop = FALSE]
stopifnot(identical(colnames(counts), names(clusters))) # para ver facilmente que los nombres sean iguales
# se crea el objeto DGEList
# primero se obtiene el grupo que hace de "meta" de las muestras
cluster <- factor(clusters[colnames(counts)], levels = 1:4, labels = c("C1","C2","C3","C4"))
# se crea el DGEList
dge <- DGEList(counts = counts, group = cluster)
# se hace un filtrado de genes con la función siguiente:
keep <- filterByExpr(dge, group = cluster)
# si se elimina un gen, recalcula los counts por muestra
dge <- dge[keep, , keep.lib.sizes = FALSE] 
dge <- calcNormFactors(dge, method = "TMM")
# se hace la matriz de diseño siguiendo el modelo cluster
design <- model.matrix(~ 0 + cluster)
colnames(design) <- levels(cluster)
# se crea el voom
v <- voom(dge, design)
# se ajusta el modelo lineal
fit <- lmFit(v, design)
# se organizan los contrastes entre clusters
contrasts <- makeContrasts(C1all = C1 - (C2+C3+C4)/3,
                           C2all = C2 - (C1+C3+C4)/3,
                           C3all = C3 - (C1+C2+C4)/3,
                           C4all = C4 - (C1+C2+C3)/3, levels = design)
# pendiente aclarar makeContrasts() y contrasts.fit()
fit_clusters <- contrasts.fit(fit, contrasts)
# treat testea cambios diferentes a un límite dado, eBayes() testea si hay diferencias respecto a 0
fit_treat <- treat(fit_clusters, lfc = 1)
genes_C1 <- topTreat(fit_treat, coef = "C1all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C2 <- topTreat(fit_treat, coef = "C2all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C3 <- topTreat(fit_treat, coef = "C3all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C4 <- topTreat(fit_treat, coef = "C4all", number = Inf, p.value = 0.05, sort.by = "logFC")
top10_C1 <- head(genes_C1, 10)
top10_C2 <- head(genes_C2, 10)
top10_C3 <- head(genes_C3, 10)
top10_C4 <- head(genes_C4, 10)

# Ahora habría que buscar mediante el entrezid o gene symbol qué genes son usando enrichGO


# 3- hacer el heatmap con esos genes 

# CUESTIONES: separar por clústers y buscar los que tienen una expresión diferencial mayor en ellos


