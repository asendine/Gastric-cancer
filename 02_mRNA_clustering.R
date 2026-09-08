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
tpm_ind <- read_excel(here("misc", "data", "tpm_filtered_log2.xlsx"))
# t(matrix) para que la disimilitud se calcule usando los genes como variables
tpm_z <- t(scale(t(tpm_ind)))


# ******************************************************************************
# clustering 1: euclidean + ward.D2
# ******************************************************************************
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")
plot(hcl, labels = FALSE, hang = -1, main = "Hierarchical clust: euc + wardd2", 
     xlab = "muestras")

# ******************************************************************************
# clustering 2: 1-Pearson + average
# ******************************************************************************
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
hcl_p <- hclust(dist_pearson, method = "average")
plot(hcl_p, labels = FALSE, hang = -1, 
     main = "Hierarchical clust: 1-Pearson + average", xlab = "muestras")

# ******************************************************************************
# clustering 3: ConsensusClusterPlus euclidean + ward.d2
# ******************************************************************************
cc_input <- as.matrix(tpm_z)
cc_results <- ConsensusClusterPlus(
  d             = cc_input,
  maxK          = 6,             # evalúa k = 2,...,6
  reps          = 500,           # 500 remuestreos
  pItem         = 0.80,          # 80% de las muestras en cada repetición
  pFeature      = 0.80,          # 80% de los genes en cada repetición
  clusterAlg    = "hc",          
  distance      = "euclidean",
  innerLinkage  = "ward.D2",     # algoritmo aplicado en cada repetición
  finalLinkage  = "average",     # agrupación final de la matriz de consenso
  seed          = 1234,
  plot          = NULL,
  writeTable    = FALSE,
  verbose       = TRUE)

# ******************************************************************************
# clustering 4: ConsensusClusterPlus + 1-pearson + average
# ******************************************************************************
cc_resultsp <- ConsensusClusterPlus(
  d             = cc_input,
  maxK          = 6,             # evalúa k = 2,...,6
  reps          = 500,           # 500 remuestreos
  pItem         = 0.80,          # 80% de las muestras en cada repetición
  pFeature      = 0.80,          # 80% de los genes en cada repetición
  clusterAlg    = "hc",          
  distance      = "pearson",
  innerLinkage  = "average",     # algoritmo aplicado en cada repetición
  finalLinkage  = "average",     # agrupación final de la matriz de consenso
  seed          = 1234,
  plot          = NULL,
  writeTable    = FALSE,
  verbose       = TRUE)

# ******************************************************************************
# ConsensusClusterPlus metrics
# ******************************************************************************
pac_results <- data.frame(
  k = 2:length(cc_results),
  PAC = sapply(2:length(cc_results), function(k) {
    mat <- cc_results[[k]]$consensusMatrix
    valores <- mat[upper.tri(mat)]
    
    mean(valores > 0.1 & valores < 0.9)
  })
)

pac_results
pac_results[which.min(pac_results$PAC), ]

icl_results <- calcICL(
  cc_results,
  plot = NULL,
  writeTable = FALSE
)

icl_results$clusterConsensus

# ******************************************************************************
pac_results <- data.frame(
  k = 2:length(cc_resultsp),
  PAC = sapply(2:length(cc_resultsp), function(k) {
    mat <- cc_resultsp[[k]]$consensusMatrix
    valores <- mat[upper.tri(mat)]
    
    mean(valores > 0.1 & valores < 0.9)
  })
)

pac_results
pac_results[which.min(pac_results$PAC), ]

icl_results <- calcICL(
  cc_resultsp,
  plot = NULL,
  writeTable = FALSE
)

icl_results$clusterConsensus


# ******************************************************************************
# NMF clustering
# ******************************************************************************
nmf_rank <- nmfEstimateRank(
  x      = tpm_ind,
  range  = 3:6,
  method = "brunet",
  nrun   = 30,
  seed   = 1234)

saveRDS(nmf_rank, file = "misc/data/nmf_rank.rds")
# ******************************************************************************
# NMF metrics
# ******************************************************************************
# cophenetic = reproductibilidad asignaciones
# dispersion = grado de consenso a valores puros 0 u 1
# evar = varianza explicada
# rss y residuals = errores respecto la matriz original
# silhouette = separación de los grupos (coef)
# sparseness = indica cada componente si esta dominado por pocos genes o muestras
nmf_rank <- readRDS("misc/data/nmf_rank.rds")
plot(nmf_rank)

consensusmap(nmf_rank$fit[["3"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["4"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["5"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["6"]], labRow = NA, labCol = NA, tracks = NA)

criteria_k <- data.frame(K=nmf_rank[["measures"]][["rank"]], 
                         Sil.coef = nmf_rank[["measures"]][["silhouette.coef"]],
                         Sil.con = nmf_rank[["measures"]][["silhouette.consensus"]])
criteria_k
