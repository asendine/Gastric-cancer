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
packages <- c()
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# libraries --------------------------------------------------------------------
libraries <- c("cBioPortalData", "TCGAbiolinks", "SummarizedExperiment", "dplyr",
               "ComplexHeatmap", "data.table", "dataframeexplorer", "devtools", 
               "cluster", "edgeR", 'limma', 'grid', 'ggplot2', 'PCAtools', 
               'ConsensusClusterPlus', 'NMF', 'readxl', 'here', 'pROC')
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
gene_info <- as.data.frame(SummarizedExperiment::rowData(rnaseq_se))
gallo_genes <- read_excel(here("misc", "otros", "2025_Gallo_claudin-low_signature", "10120_2025_1671_MOESM2_ESM.xlsx"))

# Objetivo:
# Determinar si el subtipo claudin-low definido por la firma de Gallo representa
# una clase transcriptómica estable e intrínseca de las células tumorales o un 
# fenotipo continuo condicionado por histología difusa y enriquecimiento estromal.

# Metodología artículo de mama Fougner et al.

# building the annotation table ------------------------------------------------
tpm <- assay(rnaseq_se, "tpm_unstrand")
tpm <- tpm[, substr(colnames(tpm), 14, 15) == "01", drop = FALSE]
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

summary(annotation)

# gallo_score
# seleccionamos primero solo los genes firma de Gallo (158)
# por si acaso se guarda la firma original
gallo_signature <- gallo_genes$gene_symbol
gallo_signature_old <- gallo_signature
# necesitamos filtrar los simbolos de los 158 genes de Gallo en este vector y
# aplicarlo en tpm. El problema es que Gallo tiene 18 símbolos nuevos que no constan
# en gene_info.
# se busca su símbolo antiguo usando la web https://www.genenames.org/tools/multi-symbol-checker/
# se encuentran todos excepto LOC100134259 y LOC151174.
# a través de ncbi se encuentran los dos faltantes y se añaden al archivo para que sea
# más fácil la correspondencia.
gene_check <- read.csv2(here("misc", "otros", "2025_Gallo_claudin-low_signature", "hgnc-symbol-check.csv"))
# en este archivo hay simbolo nuevo vs simbolo antiguo (el que nos interesa)
# ahora hay que seleccionar los 140 gene symbols de gallo + los 18 restantes
# se obtienen los genes a actualizar de gene_signature old
gene_check_match <- match(gallo_signature_old, gene_check$Input)
genes_to_update <- !is.na(gene_check_match)
# los genes a actualizar de gallo_signature_old se reemplazan por los antiguos
# que existen en gene_info
gallo_signature_old[genes_to_update] <- gene_check$Approved.symbol[
  gene_check_match[genes_to_update]]
# se obtienen los símboloes de los genes de las filas de tpm
tpm_gene_symbol <- gene_info$gene_name[match(rownames(tpm), gene_info$gene_id)]
# se construye el df
tpm_df <- data.frame(
  gene_id = rownames(tpm),
  gene_symbol = tpm_gene_symbol,
  tpm,
  check.names = FALSE)
# se filtra por la firma de Gallo
tpm_gallo <- tpm_df[!is.na(tpm_df$gene_symbol) & 
                         tpm_df$gene_symbol %in% gallo_signature_old, , drop = FALSE]

dim(tpm_gallo)
setdiff(gallo_signature_old, tpm_gallo$gene_symbol)
anyDuplicated(tpm_gallo$gene_symbol)

rownames(tpm_gallo) <- tpm_gallo$gene_symbol
tpm_gallo <- as.matrix(tpm_gallo_df[ , 
                                     !colnames(tpm_gallo_df) %in% 
                                       c("gene_id", "gene_symbol"), drop = FALSE])
anyNA(tpm_gallo)
anyDuplicated(rownames(tpm_gallo))
# se transforma
tpm_gallo_log <- log2(tpm_gallo + 1)
# se obtiene el z-score
tpm_gallo_z <- t(scale(t(tpm_gallo_log)))
gallo_score <- apply(tpm_gallo_z, 2, median)
# se incorpora la puntuación a annotation sabiendo que los nombres de muestra de
# tpm y rna_aliquot_id son iguales
annotation$gallo_score <- gallo_score[match(annotation$rna_aliquot_id, names(gallo_score))]
# se usa la mediana porque es lo que hacen en Londero/Gallo, la media podría ser
# una opción pero lo que buscamos es el valor típico más que el valor promedio

# validation gallo_score mayor en 56 muestras claudin low
# recogemos los valores de claudin_low_gallo dentro de annotation que NO son NA
validation_df <- annotation[!is.na(annotation$claudin_low_gallo), ]
# quedan fuera 16 muestras
aggregate(gallo_score ~ claudin_low_gallo, # cuanto difieren los grupos? resume 
          data = validation_df,            # gallo_score según grupos de claudin_low_gallo
          FUN = function(x) {              
            c(
              n = length(x),
              median = median(x),
              IQR = IQR(x)
            )
            }
          )

# son las diferencias entre los grupos claudin-low y no low significativas?
# wilcoxon para comparar distribuciones de gallo_score de los dos grupos sin 
# depender de que siga una dist normal y teniendo grupos independientes
wilcox.test(
  gallo_score ~ claudin_low_gallo,
  data = validation_df,
  exact = FALSE
)

# se confirma que el score obtenido es claramente superior en las muestras 
# claudin-low

# se calcula el AUC
roc_gallo <- pROC::roc(
  response = validation_df$claudin_low_gallo,
  predictor = validation_df$gallo_score,
  levels = c(FALSE, TRUE),
  direction = "<",
  quiet = TRUE
)

pROC::auc(roc_gallo)
# sensibilidad: claudin-low correctamente identificadas
# especificidad: no claudin-low correctamente identificadas
# se calcula cada uno en cada punto y se obtiene una métrica que resume la 
# capacidad de discriminar entre los dos grupos
# "<" indica que el gallo_score de FALSE (grupo control) debe ser menor que el
# gallo_score de TRUE (grupo claudin-low) para que la predicción sea correcta.

# genes validadores externos de la firma -> obtenidos en Gallo como genes muy 
# relacionados con claudin-low
marker_genes <- c("CLDN3", "CLDN4", "CLDN7", "CDH1", "VIM")
# se evalua la correlación entre el gallo_score y la expresión de estos genes
# se obtienen los valores de expresión de solamente los genes marcadores
marker_tpm_df <- tpm_df[tpm_df$gene_symbol %in% marker_genes, , drop = FALSE]
# se transforma a matriz quitando las columnas gene_id y gene_symbol
marker_tpm <- as.matrix(marker_tpm_df[ , !colnames(marker_tpm_df) %in% 
                                         c("gene_id", "gene_symbol"), drop = FALSE])
# se cambia el nombre de las filas a las del gene_symbol y se ordenan las filas
rownames(marker_tpm) <- marker_tpm_df$gene_symbol
marker_tpm <- marker_tpm[marker_genes, , drop = FALSE]
# se transforma
marker_log <- log2(marker_tpm + 1)
# se realiza la correlación
marker_correlations <- data.frame(gene = marker_genes, rho = NA_real_, 
                                  p_value = NA_real_)
for (i in seq_along(marker_genes)) {
  gene <- marker_genes[i]
  correlation_test <- cor.test(annotation$gallo_score, marker_log[gene, ], 
                               method = "spearman", exact = FALSE)
  marker_correlations$rho[i] <- unname(correlation_test$estimate)
  marker_correlations$p_value[i] <- correlation_test$p.value
}
marker_correlations$FDR <- p.adjust(marker_correlations$p_value, method = "BH")
marker_correlations

# se compara la expresión de estos genes marcadores entre caludin-low y no-low


# calcular stromal, immune scores y pureza

# claudin-low vs EMT no claudin-low x gallo_score

# ...?

# Breast cancer (Fougner):
# Cogen una lista de 19 genes representing only the pathognomonic gene expression
# characteristics of claudin-low tumors (manually selected on the basis of published
# characterizations of claudin-low gene expression features and advances in understanding 
# the etiological basis of claudin-low tumors).

# hacen un clustering con datos de estos genes y su expresion: Hierarchical clustering 
# using the reduced gene list was performed by complete linkage with Euclidean 
# distance as the distance metric. Para identificar el perfil claudin low poniendo 
# en comun otras variables clínicas/moleculares. We refer to tumors in this cluster
# as core claudin-low (CoreCL), while claudin-low tumors (as defined by the nine-cell
# line predictor) outside the CoreCL cluster are referred to as other claudin-low
# (OtherCL). Individual inspection of gene expression values showed that OtherCL 
# tumors displayed certain claudin-low characteristics, albeit to a lesser degree 
# than CoreCL tumors. Thus, our method for identifying claudin-low tumors primarily 
# differed from the nine-cell line predictor by filtering out a set of basal-like tumors
# with high levels of stromal and immune infiltration but without pathognomonic 
# claudin-low gene expression characteristics. 

# The significance of the core claudin-low cluster was evaluated using SigClust.

# It is therefore likely that OtherCL tumors are classified as claudin-low by the
# nine-cell line predictor due to their stromal infiltration.

# Limitaciones: The nine-cell line claudin-low predictor uses 807 genes, and Prat
# et al. acknowledge that it may inappropriately identify some tumors as claudin-low 
# solely due to stromal infiltration. We therefore considered whether a more targeted
# gene list could be used for claudin-low classification, in order to more accurately
# isolate features intrinsic to claudin-low tumors.