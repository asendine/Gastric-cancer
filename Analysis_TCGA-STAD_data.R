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
               'ConsensusClusterPlus', 'NMF')
for (i in libraries) {
  library(i, character.only = TRUE)
}

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
  file.path(tcga_dir, "Prepared", "TCGA_STAD_methylation_450K_se.rds"))
mirna_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_mirna_se.rds"))
rnaseq_se <- readRDS(file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds"))

# ******************************************************************************
# EMPEZAMOS CON RNA-SEQ*********************************************************
# ******************************************************************************

assayNames(rnaseq_se)
# unstranded: raw count
# stranded first: raw count a partir del 1a cDNA
# stranded second: raw count a partir del primer cDNA
# tpm_unstrand: Transcripts Per Million calculated using unstranded raw counts.
# fpkm_unstrand: Fragments Per Kilobase of transcript per Million mapped reads calculated from unstranded counts.
# fpkm_uq_unstrand: Upper Quartile (UQ) normalized FPKM.
tpm <- assay(rnaseq_se, "tpm_unstrand") # todas las muestras del proyecto TCGA-STAD con datos rna-seq
dim(tpm)
class(tpm)
anyNA(tpm)
# esto es 448 muestras
summary(colSums(tpm)) # paso extra cortesía de GPT, debería sumar 10^6
sum(tpm == 0, na.rm = TRUE)
# hay bastantes 0
# eliminamos las muestras con valor 0 en todos sus genes
tpm_filt <- tpm[, colSums(tpm > 0) > 0, drop = FALSE]
dim(tpm_filt)
# ahora se debería ver si todas las muestras que hay son de tumor o tejido normal
# se ve que en los nombres de las columnas, es decir muestras, la posición 14-15
# indica si es tumor (01) o tejido normal (11)
tpm_tumor <- tpm_filt[, substr(colnames(tpm_filt), 14, 15) == "01", drop = FALSE]
dim(tpm_tumor)
# hay 412 muestras de tumor, se revisa posible duplicados de muestras
patient_id <- substr(colnames(tpm_tumor), 1, 12)
sample_id  <- substr(colnames(tpm_tumor), 1, 15)
anyDuplicated(patient_id)
anyDuplicated(sample_id)
table(patient_id)[table(patient_id) > 1]
table(sample_id)[table(sample_id) > 1]
# no se elimina nada, todo parece ok
# ahora se mantendran los genes con un tpm >= 1 en el 25% de las muestras
# se calcula qué número es el 25%
n_min <- ceiling(0.25*ncol(tpm_tumor))
# se suma el valor de cada tpm por fila si es >= 1 y, se indica TRUE si el resultado es >= 103
genes_exp <- rowSums(tpm_tumor >= 1) >= n_min
# nos queda una variable con valores TRUE/FALSE por fila según si cumple o no las condiciones
tpm_filt2 <- tpm_tumor[genes_exp, , drop = FALSE]
dim(tpm_filt2)
# quedan 21k genes
# se transforman los datos a log
tpm_filt_log <- log2(tpm_filt2 + 1)

# Buscamos BACTH EFFECT
# se realiza un análisis del efecto batch simplemente por añadir calidad al workflow
# ya que sabemos por la documentación oficial que solamente hay un ligero batch effect en 
# datos de miRNA (1/4 analizados) pero no es trascendente. También se puede hacer la 
# visualización directa en la web de MDAnderson sin tener que hacerlo en R (tienen un
# web browser para ello, PCA Plus).
# A partir del nombre de las muestras se obtiene el id de placa y el TSS (origen muestra)
# que estan en las posiciones 6 y 2 respectivamente (TCGA-BR-4257-01A-01R-1131-13).
barcode_parts <- do.call(rbind, strsplit(colnames(rnaseq_se), "-", fixed = TRUE))
colData(rnaseq_se)$TSS <- factor(barcode_parts[, 2])
colData(rnaseq_se)$PlateId <- factor(barcode_parts[, 6])
batch_info <- as.data.frame(colData(rnaseq_se)[,c("TSS", "PlateId")])

pca_batch <- prcomp(t(tpm_filt_log), center = TRUE, scale. = FALSE)
variance_explained <- 100*pca_batch$sdev^2/sum(pca_batch$sdev^2)
sample_position <- match(rownames(pca_batch$x), colnames(rnaseq_se))
pca_batch_df <- data.frame(sample = rownames(pca_batch$x),
                           PC1 = pca_batch$x[, 1],
                           PC2 = pca_batch$x[, 2],
                           PlateId = batch_info$PlateId[sample_position],
                           TSS = batch_info$TSS[sample_position])

pca_background <- pca_batch_df[, c("PC1", "PC2")]

ggplot(pca_batch_df, aes(PC1, PC2)) + 
  geom_point(data = pca_background, aes(PC1, PC2), inherit.aes = FALSE, 
             color = "grey85", size = 0.7) + 
  geom_point(color = "#D55E00", size = 1.3, alpha = 0.9) + 
  facet_wrap(~ TSS, ncol = 5) +
  labs(title = "PCA: distribución de las muestras por TSS", 
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggplot(pca_batch_df, aes(PC1, PC2)) + 
  geom_point(data = pca_background, aes(PC1, PC2), inherit.aes = FALSE, 
             color = "grey85", size = 0.7) + 
  geom_point(color = "#0072B2", size = 1.3, alpha = 0.9) + 
  facet_wrap(~ PlateId, ncol = 5) +
  labs(title = "PCA: distribución de las muestras por ID de placa", 
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Ahora se buscan OUTLIERS
# para buscar outliers se realiza un PCA. Se usa PCAtools, a ver qué tal
pca_outliers <- PCAtools::pca(mat = tpm_filt_log, center = TRUE, 
                              scale = FALSE, removeVar = NULL)

# screeplot para ver % var explicada por cada componente
PCAtools::screeplot(pcaobj = pca_outliers, 
                    components = PCAtools::getComponents(pca_outliers, 1:20), 
                    title = "PCA de expresión: análisis de outliers")

# no está mal hacer un pairsplot entre las diferentes componentes, pero si hay
# muchas muestras y queremos ver muchas componentes, entonces no se verá bien
PCAtools::pairsplot(pcaobj = pca_outliers, 
                    components = PCAtools::getComponents(pca_outliers, 1:4),
                    triangle  = TRUE) # a partir de 4 la cosa empeora
# a partir de los plots anteriores, no se aprecian outliers

# ahora se puede buscar si hay genes que dominen alguna componente concreta
# se obtiene la info de genes
gene_info <- as.data.frame(rowData(rnaseq_se))
# se obtiene la posición de las filas de los genes que nos interesan con match
gene_index <- match(rownames(tpm_filt_log), rownames(rnaseq_se))
# se crea el dataset específico con ambas notaciones
gene_map <- data.frame(ensembl_id = rownames(tpm_filt_log),
                       gene_name = gene_info$gene_name[gene_index])
# a partir de aquí se crea una función para determinar los genes que más contribuyen
get_top_genes <- function(pca_object, pc, gene_map, n = 10) { # dado un pca, una componente, la leyenda y el nº de genes a comprobar
  loading_values <- pca_object$loadings[, pc] 
  # se obtienen los pesos
  top_idx <- order(abs(loading_values), decreasing = TRUE)[seq_len(n)]
  # se obtienen los n pesos por valor absoluto decreciente
  ensembl_ids <- rownames(pca_object$loadings)[top_idx]
  # se obtienen los ids de los genes con esos pesos
  data.frame(PC = pc,
             gene_name = gene_map$gene_name[match(ensembl_ids, gene_map$ensembl_id)],
             ensembl_id = ensembl_ids,
             loading = loading_values[top_idx],
             contribution_pct = 100*loading_values[top_idx]^2)
  # se hace un dataframe con la pc elegida, se busca el gene name con la leyenda indicada,
  # se indica el ensemble id, el peso obtenido y las contribuciones
  }

# aquí se puede ver qué genes dominan las componentes en este caso la 1 y 2
cont_PC1 <- get_top_genes(pca_outliers, "PC1", gene_map, n = 10)
cont_PC2 <- get_top_genes(pca_outliers, "PC2", gene_map, n = 10)
cont_PC1
cont_PC2

# ------------------------------------------------------------------------------
# DESCUBRIMIENTO DE LOS CLÚSTERS # ---------------------------------------------
# ------------------------------------------------------------------------------
# Ahora se filtra por variabilidad de expresión (median absolute deviation) buscando el 25% más variable
gene_mad <- apply(tpm_filt_log, 1, mad) # esto da un vector de mad aplicado a cada fila (1) de la matrix indicada
# igual que genes_exp pero numérico en vez de lógico
top <- ceiling(0.25*nrow(tpm_filt_log)) # se obtiene el nº de genes = 25%
ind <- order(gene_mad, decreasing = TRUE)[seq_len(top)] 
# order devuelve la posición
# seq_len crea un vector del 1 al 5349
# se obtienen las posiciones (filas) del 1 al 5349

# se limita el dataset a las filas que hay en ind
tpm_ind <- tpm_filt_log[ind, , drop=FALSE]
dim(tpm_ind)
# quedan 5349 genes

# separar estas muestras en clústers
# como queremos clusterizar muestras (columnas) habría que transponer la matriz 
# para que la disimilitud se calcule usando los genes como variables

# primero se transpone la matriz, se hace el z-score y se vuelve a transponer de nuevo.
tpm_z <- t(scale(t(tpm_ind)))
sum(is.na(tpm_z))

# clustering 1: euclidean + ward.D2. Distancia + clust
d <- dist(t(tpm_z), method = "euclidean")
hcl <- hclust(d, method = "ward.D2")
plot(hcl, labels = FALSE, hang = -1, main = "Hierarchical clust: euc + wardd2", xlab = "muestras")

# clustering 2: 1-Pearson + average, mide correlaciones entre columnas (muestras)
# se calcula la correlacion y la dist
corPearson <- cor(tpm_z, method = "pearson", use = "pairwise.complete.obs")
dist_pearson <- as.dist(1 - corPearson)
# se hace el cluster
hcl_p <- hclust(dist_pearson, method = "average")
plot(hcl_p, labels = FALSE, hang = -1, main = "Hierarchical clust: 1-Pearson + average", xlab = "muestras")

# clustering 3: ConsensusClusterPlus + euclidean/ward.d2
ccp_input <- as.matrix(tpm_z)
cc_dir <- file.path(getwd(), "output", "consensus_euclidean_wardD2")
cc_dir <- normalizePath(cc_dir, winslash = "/", mustWork = TRUE)
cc_results <- ConsensusClusterPlus(
  d             = ccp_input,
  maxK          = 6,             # evalúa k = 2,...,6
  reps          = 500,           # 500 remuestreos
  pItem         = 0.80,          # 80% de las muestras en cada repetición
  pFeature      = 0.80,          # 80% de los genes en cada repetición
  clusterAlg    = "hc",          # clustering jerárquico
  distance      = "euclidean",
  innerLinkage  = "ward.D2",     # algoritmo aplicado en cada repetición
  finalLinkage  = "average",     # agrupación final de la matriz de consenso
  seed          = 1234,
  title         = cc_dir,
  plot          = "png",
  writeTable    = FALSE,
  verbose       = TRUE
)
# clustering 4: ConsensusClusterPlus + 1-pearson + average
cc_dir2 <- file.path(getwd(), "output", "consensus_pearson_average")
cc_dir2 <- normalizePath(cc_dir2, winslash = "/", mustWork = TRUE)
cc_resultsp <- ConsensusClusterPlus(
  d             = ccp_input,
  maxK          = 6,             # evalúa k = 2,...,6
  reps          = 500,           # 500 remuestreos
  pItem         = 0.80,          # 80% de las muestras en cada repetición
  pFeature      = 0.80,          # 80% de los genes en cada repetición
  clusterAlg    = "hc",          # clustering jerárquico
  distance      = "pearson",
  innerLinkage  = "average",     # algoritmo aplicado en cada repetición
  finalLinkage  = "average",     # agrupación final de la matriz de consenso
  seed          = 1234,
  title         = cc_dir2,
  plot          = "png",
  writeTable    = FALSE,
  verbose       = TRUE
)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Al final tras probar diferentes métodos de clusterización con tpm_z:
# 1- cada método obtiene unas agrupaciones concretas
# 2- no hay una solidez en cuanto a las agrupaciones, algunos clústers son difusos
# 3- el consensus con pearson obtiene un mejor perfil con k = 3

# por ende, pienso que el problema de esta baja concordancia viene dada por el número
# de genes, el cual creo que puede ser excesivamente alto y dificulte la agrupación
# en clústers.

# cambio de planes, se prueba el método del artículo pero con las 412 muestras
# primero creamos la matrix NMF (Non-negative matrix factorization), se vuelve 
# al punto donde se creaba gene_mad y se modifica:

ind_nmf <- order(gene_mad, decreasing = TRUE)[seq_len(1500)] 
# se obtienen las posiciones (filas) del 1 al 1500
# se limita el dataset a las filas que hay en ind
tpm_nmf <- tpm_filt_log[ind_nmf, , drop=FALSE]
dim(tpm_nmf)
# quedan 1500 genes
# porqué 1500? porque si, es arbitrario y aprox lo que escogen en el estudio de 
# referencia
# evaluación del número de clústeres x algoritmo de brunet
nmf_rank <- nmfEstimateRank(
  x      = tpm_nmf,
  range  = 3:6,
  method = "brunet",
  nrun   = 30,
  seed   = 1234)

plot(nmf_rank)
nmf_rank$measures

# A partir de los resultados obtenidos lo más importante es:
# cophenetic = reproductibilidad asignaciones
# dispersion = grado de consenso a valores puros 0 u 1
# evar = varianza explicada
# rss y residuals = errores respecto la matriz original
# silhouette = separación de los grupos (coef)
# sparseness = indica cada componente si esta dominado por pocos genes o muestras

consensusmap(nmf_rank$fit[["3"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["4"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["5"]], labRow = NA, labCol = NA, tracks = NA)
consensusmap(nmf_rank$fit[["6"]], labRow = NA, labCol = NA, tracks = NA)

criteria_k <- data.frame(K=nmf_rank[["measures"]][["rank"]], 
                         Sil.coef = nmf_rank[["measures"]][["silhouette.coef"]],
                         Sil.con = nmf_rank[["measures"]][["silhouette.consensus"]])
criteria_k

# viendo las imágenes de los diferentes consensusmaps y gráficos y por las métricas obtenidas
# se podría escoger un valor de k = 3 o 6, ya que los valores diana disminuyen en k = 4 y 5.
# tampoco es que sean valores malos... es difícil, tras hacer diferentes test planteados por GPT,
# me queda la opción de hacer una consulta a claude para comparar o bien elegir
# una clusterización según mi criterio científico ya que el criterio numérico/objetivo
# es muy difuso y no es comparable a los otros métodos de clusterización realizados.

# creo que me decanto por hacer k = 6 porque a malas permite hacer una disección
# más detallada, con mayor resolución. Al final dado que el código es más o menos
# idéntico en cualquier escenario, siempre se puede cambiar y revisar los resultados.

# ------------------------------------------------------------------------------
# CARACTERIZACIÓN DE LOS CLÚSTERS # --------------------------------------------
# ------------------------------------------------------------------------------
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

# se obtienen los datos de k = 6
fit_k <- nmf_rank$fit[["6"]]
# asignación final basada en la matriz de consenso, se obtienen los clusters
clusters_k <- predict(fit_k, what = "consensus")

sil_k <- silhouette(fit_k, what = "consensus")
# se convierte sil_k en data.frame
sil_table <- data.frame(
  sample = rownames(sil_k),
  cluster = factor(sil_k[, "cluster"]),
  neighboring_cluster = factor(sil_k[, "neighbor"]),
  silhouette_width = as.numeric(sil_k[, "sil_width"]),
  row.names = NULL
)

# hora de escoger las muestras core. Criterio?




# buscar los genes diferencialmente más expresados de cada clúster (~10)
# ahora hay que hacer el análisis de expresión diferencial por cluster
# obtenemos los datos de conteos crudos para poder usar mejor edger y limma
counts <- assay(rnaseq_se, "unstranded")
dim(counts)
# se obtienen las muestras de tumor primario
counts <- counts[, names(clusters), drop = FALSE]
dim(counts)
stopifnot(identical(colnames(counts), names(clusters))) # para ver fácilmente que los nombres sean iguales

# se crea el objeto DGEList: ojo, en vez de hacer genes x muestras, se realiza el análisis
# de expresión genes x clústers
# primero se obtiene el vector categórico que hace de "meta" de las muestras
# haces factor() de clusters cogiendo la estructura de counts, indicando los niveles y etiquetas
cluster <- factor(clusters[colnames(counts)], levels = 1:4, labels = c("C1","C2","C3","C4"))
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
contrasts <- makeContrasts(C1all = C1 - (C2+C3+C4)/3,
                           C2all = C2 - (C1+C3+C4)/3,
                           C3all = C3 - (C1+C2+C4)/3,
                           C4all = C4 - (C1+C2+C3)/3, levels = design)
# contrasts.fit() coge un ajuste de un modelo lineal (fit) y una matriz con filas
# que equivalen a las columnas de los coeficientes de fit y columnas que equivalen a los contrastes?
# realiza el ajuste con contrastes
fit_clusters <- contrasts.fit(fit, contrasts)
# treat testea cambios diferentes a un límite dado (lfc > x), eBayes() testea si hay diferencias respecto a 0
fit_ebayes <- eBayes(fit_clusters)
# se obtienen los top-ranked genes más expresados diferencialmente de cada contraste
# lfc = 1 ya filtra por valor absoluto, "logFC" también
genes_C1 <- topTable(fit_ebayes, coef = "C1all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C2 <- topTable(fit_ebayes, coef = "C2all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C3 <- topTable(fit_ebayes, coef = "C3all", number = Inf, p.value = 0.05, sort.by = "logFC")
genes_C4 <- topTable(fit_ebayes, coef = "C4all", number = Inf, p.value = 0.05, sort.by = "logFC")
info_C1 <- head(genes_C1, 10) # si lfc = 0 salen muchos genes, pero si lfc = 1 solo se filtran 5!!
info_C2 <- head(genes_C2, 10) # finalmente indico 0.95 para que salgan 11 genes en C1
info_C3 <- head(genes_C3, 10) # con esto tenemos los genes en formato gene_id
info_C4 <- head(genes_C4, 10)

# 3- hacer el heatmap
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
# nos quedamos solo con las columnas que queremos y lo unimos todo
columns_keep <- c("gene_id", "gene_symbol", "gene_cluster", "logFC", "adj.P.Val")
info_C1 <- info_C1[, columns_keep]
info_C2 <- info_C2[, columns_keep]
info_C3 <- info_C3[, columns_keep]
info_C4 <- info_C4[, columns_keep]
genes_heatmap_info <- rbind(info_C1, info_C2, info_C3, info_C4)
dim(genes_heatmap_info)

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
gene_cluster <- factor(genes_heatmap_info$gene_cluster, levels = c("C1", "C2", "C3", "C4"))
sample_cluster <- factor(cluster[colnames(mat_z)], levels = c("C1", "C2", "C3", "C4"))

# se buscan duplicados y missmatches
anyDuplicated(genes_heatmap_info$gene_id)
# hay un duplicado, pero no debería ser problemático
stopifnot(identical(rownames(mat_z), genes_heatmap_info$gene_id))

# heatmap
# modificado para obtener los clusters ordenados y añadidos debajo como bloques
# por colores
cluster_colors <- c("C1" = "#E64B35","C2" = "#4DBBD5","C3" = "#00A087","C4" = "#F39B7F")

# anotación inferior 1: clústers
annotation_clusters <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors),
    labels = c("C1", "C2", "C3", "C4"),
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
subtype_colors <- c("CIN" = "#3C5488", "EBV" = "#00A087", "GS" = "#FFD700", "MSI" = "#E64B35")

# anotación inferior 2: Molecular subtype
annotation_ms <- HeatmapAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[c("C1", "C2", "C3", "C4")]),
    labels = c("C1", "C2", "C3", "C4"),
    labels_gp = gpar(col = "white",fontface = "bold")),
  Molecular_subtype = molecular_subtype,
  col = list(Molecular_subtype = subtype_colors),
  na_col = "#D9D9D9",
  annotation_height = unit(c(6, 4),"mm"),
  show_annotation_name = c(Cluster = FALSE, Molecular_subtype = TRUE),
  which = "column")

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

# Posteriormente se realizará un análisis de significación biológica para caracterizar cada clúster.
# la intención es describir cada gen dentro de cada clúster, encontrar posibles incongruencias.
# posteriormente buscar biomarcadores e hacer una matriz de muestras x biomarcador (y clúster) y plotearlo?
