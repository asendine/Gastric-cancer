# Install packages CBIO --------------------------------------------------------
bioc_packages <- c(
  "cBioPortalData",
  "TCGAbiolinks",
  "SummarizedExperiment",
  "PCAtools",
  "ConsensusClusterPlus"
)
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# install-packages-CRAN --------------------------------------------------------
packages <- c("data.table", "dataframeexplorer", "devtools", "NMF", 'writexl')
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'writexl')
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
# RNA-SEQ
# ******************************************************************************
assayNames(rnaseq_se)
# unstranded: raw count
# stranded first: raw count a partir del 1a cDNA
# stranded second: raw count a partir del primer cDNA
# tpm_unstrand: Transcripts Per Million calculated using unstranded raw counts.
# fpkm_unstrand: Fragments Per Kilobase of transcript per Million mapped reads calculated from unstranded counts.
# fpkm_uq_unstrand: Upper Quartile (UQ) normalized FPKM.

# ******************************************************************************
# Limpieza de datos y filtrado de genes
# ******************************************************************************
tpm <- assay(rnaseq_se, "tpm_unstrand")
dim(tpm)
anyNA(tpm)
# paso extra cortesía de GPT, debería sumar 10^6
summary(colSums(tpm)) 
sum(tpm == 0, na.rm = TRUE)
# eliminamos las muestras con valor 0 en todos sus genes
tpm_filt <- tpm[, colSums(tpm > 0) > 0, drop = FALSE]
dim(tpm_filt)
# muestras de tumor o tejido normal?
# tumor (01:09) o tejido normal (10:19)
tpm_tumor <- tpm_filt[, substr(colnames(tpm_filt), 14, 15) %in% sprintf("%02d", 1:9), drop = FALSE]
dim(tpm_tumor)
# duplicados de muestras
patient_id <- substr(colnames(tpm_tumor), 1, 12)
sample_id  <- substr(colnames(tpm_tumor), 1, 15)
stopifnot(anyDuplicated(patient_id) == 0)
stopifnot(anyDuplicated(sample_id) == 0)
# genes tpm >= 1 en el 25% de las muestras
n_min <- ceiling(0.25*ncol(tpm_tumor)) # nº minimo muestras
genes_exp <- rowSums(tpm_tumor >= 1) >= n_min # genes con tpm >=1 en el min de muestras
tpm_filt2 <- tpm_tumor[genes_exp, , drop = FALSE]
dim(tpm_filt2)
tpm_filt_log <- log2(tpm_filt2 + 1)
# median absolute deviation (mad) POST log
gene_mad <- apply(tpm_filt_log, 1, mad)
ind_nmf <- order(gene_mad, decreasing = TRUE)[seq_len(1500)]
tpm_ind <- tpm_filt_log[ind_nmf, , drop=FALSE]
dim(tpm_ind)
class(tpm_ind)

# se exporta el TPM filtrado y log2 a un archivo .csv
write_xlsx(as.data.frame(tpm_ind), path = "misc/data/tpm_filtered_log2.xlsx")