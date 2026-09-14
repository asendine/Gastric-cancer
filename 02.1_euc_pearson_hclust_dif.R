# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'writexl', 'here', 'fpc','NbClust')
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
# CLUSTERING
# ******************************************************************************
# data retrieving from 01_mRNA_data.R
tpm_ind <- readRDS(here("misc", "data", "tpm_ind.rds"))
# t(matrix) para que la disimilitud se calcule usando los genes como variables
tpm_z <- t(scale(t(tpm_ind)))


# ******************************************************************************
# clustering 1: euclidean + ward.D2
# ******************************************************************************
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")

# ******************************************************************************
# clustering 2: 1-Pearson + average
# ******************************************************************************
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
hcl_p <- hclust(dist_pearson, method = "average")

# ******************************************************************************
# corte - asociación de muestras a clústers
# ******************************************************************************

k <- 3

cluster_hcl <- factor(cutree(hcl, k), levels = 1:3)
cluster_hcl_p <- factor(cutree(hcl_p, k), levels = 1:3)

# ******************************************************************************
# genes diferencialmente expresados entre clústers HCL
# ******************************************************************************
counts <- assay(rnaseq_se, "unstranded")
counts <- counts[, names(cluster_hcl), drop = FALSE]
stopifnot(identical(colnames(counts), names(cluster_hcl)))
cluster <- factor(cluster_hcl[colnames(counts)], levels = 1:k, labels = paste0("C", 1:k))
dge <- DGEList(counts = counts, group = cluster)
keep <- filterByExpr(dge, group = cluster)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ 0 + cluster)
colnames(design) <- levels(cluster)

v <- voom(dge, design)
fit <- lmFit(v, design)

# ******************************************************************************
# contrastes
# ******************************************************************************

# ******************************************************************************
# heatmap
# ******************************************************************************







# ******************************************************************************
# genes diferencialmente expresados entre clústers HCL_P
# ******************************************************************************