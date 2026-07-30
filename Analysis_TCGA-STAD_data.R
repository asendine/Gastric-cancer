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

# libraries --------------------------------------------------------------------
library(cBioPortalData)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(pheatmap)

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

cbiogdc_mrna_fpkm <- read.delim(file.path(gdc_stad_dir, "data_mrna_seq_fpkm.txt"),
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

# Retreiving TCGA GDC data to classify into the different molecular subtypes


# OBJETIVO 1: CLASIFICAR LOS 443 CASOS DE GDC FROM LOS 295 DE PUB
# 1 - MEDIANTE RÉPLICA ALGORITMO (DIFÍCIL?)
# 1.1 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING Platform-Specific
# 1.1.1 - (S2)SCNA (somatic copy number alterations): SCNAs vs localización cromosoma.
        # Luego cluster jerárqico mediante GISTIC2.0 https://github.com/broadinstitute/gistic2
# 1.1.2 - ... hasta (S7)
# 1.2 - DIFINIR QUÉ VARIABLES Y MÉTODO USAN PARA EL CLUSTERING molecular data with iCluster+ (paquete bioc)

# ------------------------------------------------------------------------------

assayNames(rnaseq_se)
fpkm <- assay(rnaseq_se, "fpkm_unstrand")

prop_baja_expresion <- rowMeans(fpkm <= 0.2)
genes_keep <- prop_baja_expresion < 0.75
fpkm_filt <- fpkm[genes_keep, , drop = FALSE]

# batch effects con BLISS

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
log_fpkm_top25 <- log2(fpkm_top25 + 1)
pheatmap(
  log_fpkm_top25,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "TCGA-STAD: 25% de genes con mayor CV"
)
