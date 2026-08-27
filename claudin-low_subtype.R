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
packages <- c("data.table", "dataframeexplorer", "devtools", "NMF")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'readxl', 'here')
for (i in libraries) {
  library(i, character.only = TRUE)
}

# directory and data definition -------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")
cbio_dir <- Sys.getenv("CBIO_DATA")
pub_stad_dir <- file.path(cbio_dir, "stad_tcga_pub")
cbiopub_clin_sample <- read.delim(file.path(pub_stad_dir, "data_clinical_sample.txt"),
                                  comment.char = "#",check.names = FALSE)
cbiopub_clin_pat <- read.delim(file.path(pub_stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)
clinical_tcga <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_clinical.rds"))
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))
gallo <- read_excel(here("misc", "otros", "2025_Gallo_claudin-low_signature", "10120_2025_1671_MOESM3_ESM.xlsx"))
gallo <- gallo %>%
  mutate(MOLECULAR_SUBTYPE = recode
         (MOLECULAR_SUBTYPE, "claudin_low" = "EMT"))
gallo_final <- read_excel(here("misc", "otros", "2025_Gallo_claudin-low_signature", "10120_2025_1671_MOESM3_ESM.xlsx"))

# building the annotation table ------------------------------------------------
tpm <- assay(rnaseq_se, "tpm_unstrand")
tpm <- tpm[, substr(colnames(tpm), 14, 15) %in% sprintf("%02d", 1:9), drop = FALSE]
dim(tpm)
class(tpm)
anyNA(tpm)

gallo_match <- match(
  substr(colnames(tpm), 1, 15), substr((gallo$rna_aliquot_id), 1, 15))

gallo_final_match <- match(
  substr(colnames(tpm), 1, 15),substr((gallo_final$rna_aliquot_id), 1, 15))

annotation <- data.frame(
  rna_aliquot_id = colnames(tpm),
  sample_id = substr(colnames(tpm), 1, 15),
  patient_id = substr(colnames(tpm), 1, 12),
  tcga_subtype_original = cbiopub_clin_sample$MOLECULAR_SUBTYPE[
    match(substr(colnames(tpm), 1, 15), cbiopub_clin_sample$SAMPLE_ID)],
  acrg_subtype_deepcc = gallo$MOLECULAR_SUBTYPE[gallo_match],
  acrg_parent_subtype = gallo$MOLECULAR_SUBTYPE[gallo_match],
  claudin_low_gallo = gallo_final$MOLECULAR_SUBTYPE[gallo_final_match] == "claudin_low",
  gallo_final_subtype = gallo_final$MOLECULAR_SUBTYPE[gallo_final_match],
  histology_lauren = cbiopub_clin_pat$LAUREN_CLASS[
    match(substr(colnames(tpm), 1, 12), cbiopub_clin_pat$PATIENT_ID)],
  histology_who = cbiopub_clin_pat$WHO_CLASS[match(
    substr(colnames(tpm), 1, 12), cbiopub_clin_pat$PATIENT_ID)])

dim(annotation)

table(annotation$tcga_subtype_original, useNA = "ifany")
table(annotation$acrg_subtype_deepcc, useNA = "ifany")
table(annotation$gallo_final_subtype, useNA = "ifany")
table(annotation$claudin_low_gallo, useNA = "ifany")
table(annotation$histology_lauren, useNA = "ifany")
table(annotation$histology_who, useNA = "ifany")

anyDuplicated(annotation$sample_id)
identical(annotation$rna_aliquot_id, colnames(tpm))

annotation <- annotation %>%
  mutate(
    histology_lauren = na_if(trimws(histology_lauren), ""),
    histology_who    = na_if(trimws(histology_who), ""))
