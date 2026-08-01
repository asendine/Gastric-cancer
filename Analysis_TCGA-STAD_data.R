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

# Ahora ya se puede hacer un filtraje de los genes más representativos
# eliminamos los genes con valor 0 en todas sus muestras
tpm_tumor <- tpm_tumor[, colSums(tpm_tumor > 0) > 0, drop = FALSE]
dim(tpm_tumor)
# no hay genes con valor 0 en todas las muestras
# ahora se mantendran los genes con un tpm >= 1 en el 25% de las muestras
# se calcula qué número es el 25%
min <- ceiling(0.25*ncol(tpm_tumor))
# ahora buscamos el dataset en el que los valores sean >= 1 en 42 columnas
genes_exp <- rowSums(tpm_tumor >= 1) >= min
# nos queda un vector de valores TRUE/FALSE por fila según si cumple o no las condiciones
tpm_filt <- tpm_tumor[genes_exp, , drop = FALSE]
dim(tpm_filt)
# quedan 21k genes
# aquí se deberían evaluar si hay outliers o batch effect
# se transforman los datos a log
tpm_filt_log <- log2(tpm_filt + 1)
# Ahora se filtra por variabilidad de expresión siendo el cv = std dev/media y 
# buscando el 25% más variable
gene_mad <- apply(tpm_filt_log, 1, mad) # esto da un vector en el hay las Median
# Absolute Deviation de cada fila (gen) (1) de la matrix indicada
top <- ceiling(0.25*nrow(tpm_filt_log)) # se obtiene el nº de genes = 25%
ind <- order(gene_mad, decreasing = TRUE)[seq_len(top)] 
# order devuelve la posición
# seq_len crea un vector del 1 al 5349
# se obtienen las posiciones (filas) del 1 al 5349

# se limita el dataset a las filas que hay en ind
tpm_filt <- tpm_filt_log[ind, , drop=FALSE]
dim(tpm_filt)
# quedan 5349 genes

library(devtools)
install_github("jokergoo/ComplexHeatmap")



pheatmap(
  tpm_filt,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "TCGA-STAD: 25% de genes con mayor MAD"
)


media_fpkm <- rowMeans(fpkm_filt, na.rm = TRUE)
keep_mean <- is.finite(media_fpkm) & media_fpkm >= 10
fpkm_mean10 <- fpkm_filt[keep_mean, , drop = FALSE]
media_mean10 <- media_fpkm[keep_mean]
sd_fpkm <- apply(fpkm_mean10, 1, sd, na.rm = TRUE)
cv_fpkm <- sd_fpkm / media_mean10
keep_cv <- is.finite(cv_fpkm)
fpkm_cv <- fpkm_mean10[keep_cv, , drop = FALSE]
cv_fpkm <- cv_fpkm[keep_cv]
orden_cv <- order(cv_fpkm, decreasing = TRUE)
n_top <- ceiling(0.25 * nrow(fpkm_cv))
top_idx <- orden_cv[seq_len(n_top)]
fpkm_top25 <- fpkm_cv[top_idx, , drop = FALSE]

