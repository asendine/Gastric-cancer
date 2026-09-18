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

k_values <- 2:10

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
  maxK          = 10,             # evalúa k = 2,...,6
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

ml <- vector("list", length(cc_results))

for (k in 2:length(cc_results)) {
  ml[[k]] <- cc_results[[k]]$consensusMatrix
}

ConsensusClusterPlus:::CDF(ml)

# ******************************************************************************
# clustering 4: ConsensusClusterPlus + 1-pearson + average
# ******************************************************************************
cc_resultsp <- ConsensusClusterPlus(
  d             = cc_input,
  maxK          = 10,             # evalúa k = 2,...,6
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

mlp <- vector("list", length(cc_resultsp))

for (k in 2:length(cc_resultsp)) {
  mlp[[k]] <- cc_resultsp[[k]]$consensusMatrix
}

ConsensusClusterPlus:::CDF(mlp)

# ******************************************************************************
# ConsensusClusterPlus metrics
# ******************************************************************************
# PAC
# ICL
# CDF (DELTA AREA)

pac_results <- data.frame(
  k = 2:length(cc_results),
  PAC = sapply(2:length(cc_results), function(k) {
    mat <- cc_results[[k]]$consensusMatrix # obtiene la matriz de consenso
    valores <- mat[upper.tri(mat)] # obtiene los valores de la diagonal principal
    mean(valores > 0.1 & valores < 0.9) # calcula la media de los valores ambiguos
  }
  )
)
pac_resultsp <- data.frame(
  k = 2:length(cc_resultsp),
  PAC = sapply(2:length(cc_resultsp), function(k) {
    mat <- cc_resultsp[[k]]$consensusMatrix
    valores <- mat[upper.tri(mat)]
    mean(valores > 0.1 & valores < 0.9)
  }
  )
)

# ******************************************************************************
# GRÁFICO PAC
# ******************************************************************************
pac_plot <- rbind(transform(pac_results,  method = "Euclidean"),
                  transform(pac_resultsp, method = "Pearson"))

ggplot(pac_plot, aes(x = k, y = PAC, group = method, shape = method)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = sort(unique(pac_plot$k))) +
  labs(x = "k", y = "PAC", shape = NULL) +
  theme_classic()

# ******************************************************************************
# GRÁFICO ICL
# ******************************************************************************
k_selected <- c(3, 4, 5)

icl_euclidean <- as.data.frame(icl_results$clusterConsensus)
icl_pearson   <- as.data.frame(icl_resultsp$clusterConsensus)
icl_euclidean <- icl_euclidean[icl_euclidean$k %in% k_selected,]
icl_pearson <- icl_pearson[icl_pearson$k %in% k_selected,]
icl_euclidean$method <- "Euclidean + Ward.D2"
icl_pearson$method   <- "Pearson + Average"

icl_plot <- rbind(icl_euclidean,icl_pearson)

ggplot(icl_plot, aes(x = factor(k), y = clusterConsensus, colour = method, group = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.35)) +
  geom_text(aes(label = cluster), position = position_dodge(width = 0.35), 
            hjust = -1, size = 4, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "ICL values per k and method", x = "k", y = "Cluster consensus", 
       colour = "Method") +
  theme_classic()

# ******************************************************************************
# NMF clustering
# ******************************************************************************
nmf_rank <- nmfEstimateRank(
  x      = tpm_ind,
  range  = 3:8,
  method = "brunet",
  nrun   = 30,
  seed   = 1234)

saveRDS(nmf_rank, file = "misc/data/nmf_rank.rds")

# ******************************************************************************
# NMF metrics
# ******************************************************************************
# cophenetic: measures how reliably the same samples are assigned to the same cluster.
# sil.coef: evalúa la separación de las muestras en la matriz de coeficientes.
# sil. con: silueta del consenso (1 − consenso). Compara la proximidad de una 
# muestra a su propio clúster con la proximidad al clúster alternativo más cercano
nmf_rank <- readRDS("misc/data/nmf_rank.rds")
plot(nmf_rank)

consensusmap(nmf_rank$fit[["4"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["6"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["7"]], labRow = NA, labCol = NA, tracks = NA)

criteria_k <- data.frame(K=nmf_rank[["measures"]][["rank"]], 
                         Sil.coef = nmf_rank[["measures"]][["silhouette.coef"]],
                         Sil.con = nmf_rank[["measures"]][["silhouette.consensus"]])
criteria_k

# ******************************************************************************
# GRÁFICO SILUETAS NMF
# ******************************************************************************
k_selected <- c(4, 6, 7)

sil_plot <- do.call(rbind, lapply(k_selected, function(k) {
  sil <- silhouette(nmf_rank$fit[[as.character(k)]], what = "consensus")
  df <- as.data.frame(sil[, c("cluster", "sil_width")])
  df <- df[order(df$cluster, -df$sil_width), ]
  df$position <- seq_len(nrow(df))
  df$panel <- sprintf("k = %s | Silueta media = %.3f", k, mean(df$sil_width))
  df$mean_sil <- mean(df$sil_width)
  df
}))

sil_plot$panel <- factor(sil_plot$panel, levels = unique(sil_plot$panel))


ggplot(sil_plot, aes(x = position, y = sil_width, fill = factor(cluster))) +
  geom_col(width = 1) +
  geom_hline(data = unique(sil_plot[, c("panel", "mean_sil")]), 
             aes(yintercept = mean_sil), colour = "red", linetype = "dashed") +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = "Muestras", y = "Anchura de silueta", fill = "Clúster") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

# ******************************************************************************
# número de clústers
# ******************************************************************************

# hcl
# pearson -> k = 3 o 4 aprox igual

# ConsensusClusterPlus
# pearson -> k = 4

# NMF clustering
# k = 4 o 6

