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
# CARACTERIZACIÓN DE LOS CLÚSTERS
# ******************************************************************************
# el siguiente paso es encontrar las "core samples", aquellas con anchura de silueta
# más positiva ~mayor similaridad con su cluster.
# "Samples most representative of the clusters, hereby called core samples were 
# identified based on positive silhouette width, indicating higher similarity to 
# their own class than to any other class member. Core samples were used to select 
# differentially expressed marker genes for each subtype by comparing the subclass 
# versus the other subclasses, using Student's t-test." La realidad es que usando 
# cualquier valor de k se han obtenido valores extremadamente buenos de anchura de
# silueta. Aquí surge otro dilema, plantear un criterio arbitrario para filtrar muestras
# o dejarlas todas y ya está... no me parece mal en un contexto exploratorio
# obtener solamente las muestras que más definan un cluster sin pasarnos. El criterio
# silueta > 0 no discrimina nada, empezaría a ser interesante con sil > 0.6, en el que
# se eliminarían algunas muestras, quizá un 30%, que ya está bien.

nmf_rank <- readRDS("misc/data/nmf_rank.rds")

# se obtienen los datos de k = 6
fit <- nmf_rank$fit[["6"]]
# asignación final basada en la matriz de consenso, se obtienen los clusters
clusters <- predict(fit, what = "consensus")

names(clusters)  # barcodes
clusters         # números de clúster

sil <- silhouette(fit, what = "consensus")
# se convierte sil_k en data.frame
sil_table <- data.frame(
  sample = rownames(sil),
  cluster = factor(sil[, "cluster"]),
  neighboring_cluster = factor(sil[, "neighbor"]),
  silhouette_width = as.numeric(sil[, "sil_width"]),
  row.names = NULL
)

counts <- assay(rnaseq_se, "unstranded")
dim(counts)
# se obtienen las muestras de tumor primario
counts <- counts[, names(clusters), drop = FALSE]
dim(counts)
stopifnot(identical(colnames(counts), names(clusters))) # para ver fácilmente que los nombres sean iguales
# se crea el objeto DGEList: ojo, en vez de hacer genes x muestras, se realiza el análisis
# de expresión genes x clústers
# primero se obtiene el vector categórico que hace de "meta" de las muestras
# se hace factor() de clusters cogiendo la estructura de counts, indicando los niveles y etiquetas
cluster <- factor(clusters[colnames(counts)], levels = 1:6, labels = paste0("C", 1:6))
dge <- DGEList(counts = counts, group = cluster)
# ---------
# se hace un filtrado de genes preventivo con filterByExpr(). Filtrado por defecto, sin más.
keep <- filterByExpr(dge, group = cluster)
# si se elimina un gen, recalcula los tamaños de biblioteca como la suma de los genes conservados
# esto es cosa de GPT, nunca lo he usado:
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
# ---------
# se hace la matriz de diseño siguiendo el modelo cluster
design <- model.matrix(~ 0 + cluster)
head(design)
# se cambian los nombres de las columnas
colnames(design) <- levels(cluster)
head(design)
# se crea el voom para hacer el modelo lineal (dge + matriz coeficientes)
# v contiene logCPM a partir de raw counts
v <- voom(dge, design)
# se ajusta el modelo lineal
fit <- lmFit(v, design)
# se organizan los contrastes entre clusters, se compara la expresión de cada clúster
# con la media de expresiones del resto... es una forma pero hay otras
# sintaxis es: diferentes contrastes como un vector + los niveles = matriz de diseño.
contrasts_mat <- makeContrasts(C1all = C1 - (C2+C3+C4+C5+C6)/5,
                           C2all = C2 - (C1+C3+C4+C5+C6)/5,
                           C3all = C3 - (C1+C2+C4+C5+C6)/5,
                           C4all = C4 - (C1+C2+C3+C5+C6)/5,
                           C5all = C5 - (C1+C2+C3+C4+C6)/5,
                           C6all = C6 - (C1+C2+C3+C4+C5)/5, levels = design)
# contrasts.fit() coge un ajuste de un modelo lineal (fit) y una matriz con filas
# que equivalen a las columnas de los coeficientes de fit y columnas que equivalen a los contrastes?
# realiza el ajuste con contrastes
fit_clusters <- contrasts.fit(fit, contrasts_mat)
# treat testea cambios diferentes a un límite dado (lfc > x), eBayes() testea si hay diferencias respecto a 0
fit_ebayes <- eBayes(fit_clusters)
# se obtienen los top-ranked genes más expresados diferencialmente de cada contraste
# lfc = 1 ya filtra por valor absoluto, "logFC" también
genes_C1 <- topTable(fit_ebayes, coef = "C1all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C2 <- topTable(fit_ebayes, coef = "C2all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C3 <- topTable(fit_ebayes, coef = "C3all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C4 <- topTable(fit_ebayes, coef = "C4all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C5 <- topTable(fit_ebayes, coef = "C5all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C6 <- topTable(fit_ebayes, coef = "C6all", number = Inf, p.value = 0.05, sort.by = "logFC")
info_C1 <- head(genes_C1, 10) # si lfc = 0 salen muchos genes, pero si lfc = 1 solo se filtran 5!!
info_C2 <- head(genes_C2, 10) # finalmente indico 0.95 para que salgan 11 genes en C1
info_C3 <- head(genes_C3, 10) # con esto tenemos los genes en formato gene_id
info_C4 <- head(genes_C4, 10)
info_C5 <- head(genes_C5, 10)
info_C6 <- head(genes_C6, 10)

# 3- hacer el heatmap
gene_info <- as.data.frame(rowData(rnaseq_se))
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
# C5 ----------------------------------------------------------
info_C4$gene_id <- rownames(info_C4)
info_C4$gene_symbol <- gene_info$gene_name[match(info_C4$gene_id, rownames(gene_info))]
info_C4$gene_cluster <- "C4"
head(info_C4)
# C4 ----------------------------------------------------------
info_C4$gene_id <- rownames(info_C4)
info_C4$gene_symbol <- gene_info$gene_name[match(info_C4$gene_id, rownames(gene_info))]
info_C4$gene_cluster <- "C4"
head(info_C4)
# C5 ----------------------------------------------------------
info_C5$gene_id <- rownames(info_C5)
info_C5$gene_symbol <- gene_info$gene_name[match(info_C5$gene_id, rownames(gene_info))]
info_C5$gene_cluster <- "C5"
head(info_C5)
# C6 ----------------------------------------------------------
info_C6$gene_id <- rownames(info_C6)
info_C6$gene_symbol <- gene_info$gene_name[match(info_C6$gene_id, rownames(gene_info))]
info_C6$gene_cluster <- "C6"
head(info_C6)

#nos quedamos con las columnas que queremos y lo unimos todo
columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")
info_C1 <- info_C1[, columns_keep]
info_C2 <- info_C2[, columns_keep]
info_C3 <- info_C3[, columns_keep]
info_C4 <- info_C4[, columns_keep]
info_C4 <- info_C4[, columns_keep]
info_C5 <- info_C5[, columns_keep]
info_C6 <- info_C6[, columns_keep]
genes_heatmap_info <- rbind(info_C1, info_C2, info_C3, info_C4, info_C5, info_C6)

# Segundo se construye la matriz de expresión
# el objeto "v" creado por voom contiene la relación entre el gene_id
# y las muestras. Los valores son logCPM. $weights contiene los pesos
gene_ids <- genes_heatmap_info$gene_id
# se filtran los 40 genes que nos interesan
mat_heatmap <- v$E[gene_ids, , drop = FALSE]
# y se calcula el z-score de cada valor
# básicamente para cada elemento en mat_heatmap se calcula la media y la sd
# entonces se calcula el z-score de cada valor con la media y la sd y se sustituye
# en mat_z:
mat_z <- mat_heatmap
for (i in 1:nrow(mat_heatmap)) {
  gene_mean <- mean(mat_heatmap[i, ])
  gene_sd <- sd(mat_heatmap[i, ])
  mat_z[i, ] <- (mat_heatmap[i, ] - gene_mean)/gene_sd
}
# por último se crea la agrupación por clusters para las filas y las columnas
gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = c("C1", "C2", "C3", "C4", "C5", "C6"))
sample_cluster <- factor(cluster[colnames(mat_z)], levels = c("C1", "C2", "C3", "C4", "C5", "C6"))

# se buscan duplicados y missmatches
anyDuplicated(genes_heatmap_info$gene_id)

# heatmap
# modificado para obtener los clusters ordenados y añadidos debajo como bloques
# por colores
cluster_colors <- c("C1" = "#FF0000","C2" = "blue","C3" = "#00FF00","C4" = "#FFFF00","C5" = "#FFA500","C6" = "#A020F0")
# anotación inferior 1: clústers
annotation_clusters <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors),
    labels = c("C1", "C2", "C3", "C4", "C5", "C6"),
    labels_gp = gpar(col = "white", fontface = "bold")),
  which = "column")

# finalmente se relaciona con los molecular subtypes obtenidos en las 290 y pico muestras
sample_id_mat <- substr(colnames(mat_z), 1, 15)
posicion <- match(sample_id_mat, cbiopub_clin_sample$SAMPLE_ID)
sum(is.na(posicion))
# hay 138 NAs, es decir 274 muestras etiquetadas con subtype de las 412, pero en 
# cbio constan 295 muestras etiquetadas. Hay 21 muestras que no estan entre las 412?
molecular_subtype <- cbiopub_clin_sample$MOLECULAR_SUBTYPE[posicion]
names(molecular_subtype) <- colnames(mat_z)
head(molecular_subtype)
subtype_colors <- c("CIN" = "darkorchid4", "EBV" = "#00A087", "GS" = "#B8860B", "MSI" = "#FF69B4")

# anotación inferior 2: Molecular subtype
annotation_ms <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[c("C1", "C2", "C3", "C4", "C5", "C6")]),
    labels = c("C1", "C2", "C3", "C4", "C5", "C6"),
    labels_gp = gpar(col = "white",fontface = "bold")),
  Molecular_subtype = molecular_subtype,
  col = list(Molecular_subtype = subtype_colors),
  na_col = "white",
  annotation_height = unit(c(6, 4),"mm"),
  show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE),
  which = "column")

# cor Pearson sobre filas
row_corPearson <- cor(t(mat_z), method = "pearson")
row_dist_pearson <- as.dist(1 - row_corPearson)
# se hace el cluster
row_hcl_p <- hclust(row_dist_pearson, method = "average")

# cluster genes (separacion genes por tipo de cluster)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  row_split = gene_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)

# se busca el vector de clusters según orden de genes
marker_cluster <- genes_heatmap_info$gene_cluster[match(rownames(mat_z), 
                                                        genes_heatmap_info$gene_id)]
marker_cluster <- factor(marker_cluster, levels = paste0("C", 1:6))
row_annotation <- rowAnnotation(Clusters = marker_cluster, 
                                col = list(Clusters = c(
      C1 = "red",
      C2 = "blue",
      C3 = "green",
      C4 = "yellow",
      C5 = "orange",
      C6 = "purple")))

# row clustering 1-pearson + average (como se expresan en general)
Heatmap(
  mat_z,
  name = "Z-score",
  row_labels = genes_heatmap_info$gene_symbol,
  row_names_gp = gpar(fontsize = 7),
  column_split = sample_cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = row_hcl_p,
  row_dend_reorder = FALSE,
  row_split = NULL,
  left_annotation = row_annotation,
  bottom_annotation = annotation_ms,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)

# pendiente estudiar el paquete NMF para entenderlo. Estudiar como funciona el metodo brunet.


# Posteriormente se realizará un análisis de significación biológica para caracterizar cada clúster.


