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
# hclust metrics
# ******************************************************************************

# silhouette (-1 a 1): distancia media de la muestra en su cluster vs otro cluster.
# Media por clúster o global.

# Índice de Dunn: La distancia más pequeña entre dos puntos que pertenecen a 
# clústeres diferentes / la distancia más grande entre dos puntos dentro del 
# mismo clúster. Mayor = mejor.

# Calinski–Harabasz: mide el grado de separación entre grupos y si estos son 
# compactos. Mayor = mejor.

k_values <- 2:6

arboles <- list(euclidean_wardD2 = hcl, pearson_average = hcl_p)

distancias <- list(euclidean_wardD2 = d, pearson_average = dist_pearson)

metricas <- data.frame()

silhouette_por_cluster <- data.frame()

for (metodo in names(arboles)) {
  for (k in k_values) {
    
    grupos <- cutree(arboles[[metodo]], k = k) # genera los k clusters
    # paquete fpc, cluster.stats obtiene las métricas que necesito
    res <- cluster.stats(d = distancias[[metodo]],
                         clustering = grupos,
                         wgap = FALSE, sepindex = FALSE) # no necesarios 
    # se obtienen las métricas
    metricas <- rbind(metricas, data.frame(metodo = metodo,
                                           k = k,
                                           silhouette = res$avg.silwidth,
                                           dunn = res$dunn,
                                           calinski_harabasz = res$ch))
    # y las siluetas por cluster
    silhouette_por_cluster <- rbind(silhouette_por_cluster, data.frame(
      metodo = metodo,
      k = k,
      cluster = seq_len(k),
      n = res$cluster.size,
      silhouette = as.numeric(res$clus.avg.silwidths)))
  }
}

metricas

silhouette_por_cluster

# ******************************************************************************
# usando NbClust (lo mismo pero más automátizado)
# ******************************************************************************

indices <- c("silhouette", "dunn", "cindex", "mcclain")

nb_euc <- list()
nb_pearson <- list()

for (indice in indices) {
  
  nb_euc[[indice]] <- NbClust(data = NULL,
                              diss = d,
                              distance = NULL,
                              min.nc = 3, max.nc = 6,
                              method = "ward.D2",
                              index = indice)
  
  nb_pearson[[indice]] <- NbClust(data = NULL,
                                  diss = dist_pearson,
                                  distance = NULL,
                                  min.nc = 3, max.nc = 6,
                                  method = "average",
                                  index = indice)
}

# Valores de los indices para cada k
metricas_nb_euc <- sapply(nb_euc, "[[", "All.index")
metricas_nb_pearson <- sapply(nb_pearson, "[[", "All.index")

# k recomendado por cada indice y valor obtenido
mejores_k_euc <- sapply(nb_euc, "[[", "Best.nc")
mejores_k_pearson <- sapply(nb_pearson, "[[", "Best.nc")


metricas_nb_euc
metricas_nb_pearson

mejores_k_euc
mejores_k_pearson

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
# PAC = Proportion of Ambiguous Clustering, es la proporción de valores de la 
# matriz de consenso que se encuentran entre 0.1 y 0.9, es decir, que no son ni 
# 0 ni 1. Un valor bajo de PAC indica que la mayoría de las muestras se asignan 
# a un solo cluster con alta probabilidad, lo que sugiere una buena estabilidad 
# del clustering.
# miramos PACs (bajo = mejor)
pac_results <- data.frame(
  k = 2:length(cc_results),
  PAC = sapply(2:length(cc_results), function(k) {
    mat <- cc_results[[k]]$consensusMatrix # obtiene la matriz de consenso
    valores <- mat[upper.tri(mat)] # obtiene los valores de la diagonal principal
    mean(valores > 0.1 & valores < 0.9) # calcula la media de los valores ambiguos
  }
  )
)

pac_results
# el icl es un índice de estabilidad de clustering que combina la información 
# de la matriz de consenso y la asignación de clusters. Un valor más alto de ICL
# indica una mayor estabilidad del clustering.
icl_results <- calcICL(cc_results, plot = NULL, writeTable = FALSE)
icl_results$clusterConsensus

# ******************************************************************************
pac_resultsp <- data.frame(
  k = 2:length(cc_resultsp),
  PAC = sapply(2:length(cc_resultsp), function(k) {
    mat <- cc_resultsp[[k]]$consensusMatrix
    valores <- mat[upper.tri(mat)]
    mean(valores > 0.1 & valores < 0.9)
  }
  )
)

pac_resultsp

pac_resultsp[which.min(pac_resultsp$PAC), ]

icl_resultsp <- calcICL(cc_resultsp, plot = NULL, writeTable = FALSE)
icl_resultsp$clusterConsensus

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

# ******************************************************************************
# número de clústers
# ******************************************************************************

# hcl y hcl_p
# euclidean -> k = 3 o quizá 5 
# pearson -> k = 3 o 4 aprox igual

# ConsensusClusterPlus
# euc -> k = 6 (0.5153072)
# pearson -> k = 6 (0.4740510)
# cuanto más alto el valor de k menor es el PAC. Esto no cuadra mucho, pdnte investigar.

# NMF clustering
# k = 3 - Sil.coef(0.6362434) - Sil.con(0.9277994)

# procederemos con k = 3 por ahora

# ******************************************************************************
# corte - asociación de muestras a clústers
# ******************************************************************************

k <- 3

