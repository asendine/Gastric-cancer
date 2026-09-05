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
               'ConsensusClusterPlus', 'NMF', 'readxl', 'here', 'pROC',
               'tidyestimate', 'tidyr')
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
gallo_parent <- gallo %>%
  mutate(MOLECULAR_SUBTYPE = recode(MOLECULAR_SUBTYPE, "claudin_low" = "EMT"))
gene_info <- as.data.frame(SummarizedExperiment::rowData(rnaseq_se))
gallo_genes <- read_excel(here("misc", "otros", "2025_Gallo_claudin-low_signature", "10120_2025_1671_MOESM2_ESM.xlsx"))

# ******************************************************************************
# Objetivo:
# Determinar si el subtipo claudin-low definido por la firma de Gallo representa
# una clase transcriptómica estable e intrínseca de las células tumorales o un 
# fenotipo continuo condicionado por histología difusa y enriquecimiento estromal.

# Metodología artículo de mama Fougner et al.
# ******************************************************************************
# building the annotation table
# ******************************************************************************
tpm <- assay(rnaseq_se, "tpm_unstrand")
tpm <- tpm[, substr(colnames(tpm), 14, 15) == "01", drop = FALSE]
dim(tpm)
class(tpm)
if (anyNA(tpm) || any(!is.finite(tpm))) {
  stop("tpm contiene valores NA, NaN o infinitos.")}
if (anyDuplicated(rownames(tpm))) {
  stop("Los gene_id de tpm no son únicos.")}
if (any(tpm < 0)) {stop("Los valores TPM no pueden ser negativos.")}


# match clasificación gallo
gallo_match <- match( 
  substr(colnames(tpm), 1, 15), substr((gallo$rna_aliquot_id), 1, 15))
# match clasificación ACRG
gallo_parent_match <- match( 
  substr(colnames(tpm), 1, 15), substr((gallo_parent$rna_aliquot_id), 1, 15))

acrg_subtype <- gallo_parent$MOLECULAR_SUBTYPE[gallo_parent_match]
gallo_subtype <- gallo$MOLECULAR_SUBTYPE[gallo_match]

annotation <- data.frame(
  rna_aliquot_id = colnames(tpm),
  sample_id = substr(colnames(tpm), 1, 15),
  patient_id = substr(colnames(tpm), 1, 12),
  tcga_subtype_original = cbiopub_clin_sample$MOLECULAR_SUBTYPE[
    match(substr(colnames(tpm), 1, 15), cbiopub_clin_sample$SAMPLE_ID)],
  acrg_parent_subtype = acrg_subtype,
  claudin_low_gallo = gallo_subtype == "claudin_low",
  gallo_subtype = gallo_subtype,
  histology_lauren = cbiopub_clin_pat$LAUREN_CLASS[
    match(substr(colnames(tpm), 1, 12), cbiopub_clin_pat$PATIENT_ID)],
  histology_who = cbiopub_clin_pat$WHO_CLASS[match(
    substr(colnames(tpm), 1, 12), cbiopub_clin_pat$PATIENT_ID)])

dim(annotation)

table(annotation$tcga_subtype_original, useNA = "ifany")
table(annotation$acrg_parent_subtype, useNA = "ifany")
table(annotation$gallo_subtype, useNA = "ifany")
table(annotation$claudin_low_gallo, useNA = "ifany")
table(annotation$histology_lauren, useNA = "ifany")
table(annotation$histology_who, useNA = "ifany")

stopifnot(anyDuplicated(annotation$sample_id) == 0L, 
          identical(annotation$rna_aliquot_id, colnames(tpm)))

annotation <- annotation %>%
  mutate(
    histology_lauren = na_if(trimws(histology_lauren), ""),
    histology_who    = na_if(trimws(histology_who), ""))

summary(annotation)

# ******************************************************************************
# gallo_score into annotation
# ******************************************************************************
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
# los genes a cambiar de gallo_signature_old (genes_to_update) se reemplazan por 
# los antiguos que existen en gene_info (gene_check$Approved.symbol)
gallo_signature_old[genes_to_update] <- gene_check$Approved.symbol[
  gene_check_match[genes_to_update]]
# se obtienen los símbolos de los genes de las filas de tpm
tpm_gene_symbol <- gene_info$gene_name[match(rownames(tpm), gene_info$gene_id)]

# se construye el df
tpm_df <- data.frame(
  gene_id = rownames(tpm),
  gene_symbol = tpm_gene_symbol,
  tpm,
  check.names = FALSE)

if (anyDuplicated(tpm_df$gene_id)) {stop("tpm_df contiene gene_id duplicados.")}

# se filtra por la firma de Gallo
tpm_gallo <- tpm_df[!is.na(tpm_df$gene_symbol) & 
                         tpm_df$gene_symbol %in% gallo_signature_old, , drop = FALSE]

dim(tpm_gallo)

rownames(tpm_gallo) <- tpm_gallo$gene_symbol
tpm_gallo <- as.matrix(tpm_gallo[ , 
                                     !colnames(tpm_gallo) %in% 
                                       c("gene_id", "gene_symbol"), drop = FALSE])
# se transforma
tpm_gallo_log <- log2(tpm_gallo + 1)
gene_sd <- apply(tpm_gallo_log, 1, sd)
if (any(!is.finite(gene_sd) | gene_sd == 0)) {
  stop("Hay genes de la firma sin variabilidad entre muestras.")}

# se obtiene el z-score
tpm_gallo_z <- t(scale(t(tpm_gallo_log)))
stopifnot(!anyNA(tpm_gallo_z), all(is.finite(tpm_gallo_z)))
# se crea la puntuación de Gallo
gallo_score <- apply(tpm_gallo_z, 2, median)
# se incorpora la puntuación a annotation sabiendo que los nombres de muestra de
# tpm y rna_aliquot_id son iguales
annotation$gallo_score <- gallo_score[match(annotation$rna_aliquot_id, names(gallo_score))]
# se usa la mediana porque es lo que hacen en Londero/Gallo, la media podría ser
# una opción pero lo que buscamos es el valor típico más que el valor promedio

# ******************************************************************************
# Checks gallo_score
# ******************************************************************************
# gallo_score en grupos claudin-low/no-low diferencias entre grupos
# usando agregate + wilcoxon
# ******************************************************************************
# se construye un df sin las muestras NA en claudin_low_gallo
validation_df <- annotation[!is.na(annotation$claudin_low_gallo), ]
# quedan fuera 16 muestras
aggregate(gallo_score ~ claudin_low_gallo, # cuanto difieren los grupos? resume 
          data = validation_df,            # gallo_score según grupos de claudin_low_gallo
          FUN = function(x) {              
            c(n = length(x),
              median = median(x),
              IQR = IQR(x))
            })

# son las diferencias entre los grupos claudin-low y no low significativas?
# wilcoxon para comparar distribuciones de gallo_score de los dos grupos sin 
# depender de que siga una dist normal y teniendo grupos independientes
wilcox.test(gallo_score ~ claudin_low_gallo, data = validation_df,  exact = FALSE)

# se confirma que el score obtenido es claramente superior en las muestras 
# claudin-low

# ******************************************************************************
# AUC: dado un gallo_score, qué tan bien se identificaría correctamente en 
# claudin-low o no-low
# ******************************************************************************
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

# ******************************************************************************
# Marker genes 1: correlation entre genes marker y gallo_score
# ******************************************************************************
# genes check de la firma -> obtenidos en Gallo como genes muy relacionados con 
# claudin-low
marker_genes <- c("CLDN3", "CLDN4", "CLDN7", "CDH1", "VIM")
# se comprueba que ninguno de los genes de la firma se encuentra entre los 
# markers
intersect(marker_genes, gallo_signature_old)
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
  # se esta evaluando el primer valor de annotation (gallo_score del gen 1) con 
  # el primer valor del gen 1 de marker_log porque las columnas siguen el mismo 
  # orden que las filas de annotation.
  correlation_test <- cor.test(annotation$gallo_score, marker_log[gene, ], 
                               method = "spearman", exact = FALSE)
  # guarda el número limpio:
  marker_correlations$rho[i] <- unname(correlation_test$estimate)
  marker_correlations$p_value[i] <- correlation_test$p.value
}
marker_correlations$FDR <- p.adjust(marker_correlations$p_value, method = "BH")
marker_correlations
# rho: es la correlacion de spearman que evalua la fuerza y dirección de la 
# relación entre dos variables.
# pearson -> relacion lineal
# spearman -> relacion monotona (no necesariamente lineal)
# kendall -> no se usa normalmente en transcriptómica
# BH -> Benjamini-Hochberg corrección que evita obtener falsos positivos (FDR) 
# entre los resultados significativos cuando estos son muchos y pueden ser por 
# azar. Bonferroni aplica una corrección mucho más ajustada y puede eliminar 
# resultados buenos.

# ******************************************************************************
# Marker genes 2: comparison between claudin-low and non-low
# comparación expresión entre genes marker vs claudin-low y no-low.
# ******************************************************************************
# tpm_log de cada gen marcador vs 56 muestras claudin-low y el tpm_log de cada 
# gen marcador con las 340 muestras no-low usando Wilcoxon.
sample_match <- match(annotation$rna_aliquot_id, colnames(marker_log))
marker_group_df <- data.frame(rna_aliquot_id = annotation$rna_aliquot_id,
                              claudin_low_gallo = annotation$claudin_low_gallo,
                              t(marker_log[, sample_match, drop = FALSE]),
                              check.names = FALSE)
# se eliminan los NAs
marker_group_df <- marker_group_df[!is.na(marker_group_df$claudin_low_gallo), 
                                   , drop = FALSE]
dim(marker_group_df)

# se construye una tabla con los resultados
marker_group_results <- data.frame(
  gene = marker_genes,
  median_non_low = NA_real_,
  IQR_non_low = NA_real_,
  median_claudin_low = NA_real_,
  IQR_claudin_low = NA_real_,
  median_difference = NA_real_,
  p_value = NA_real_)

# comparación de cada gen vs fenotipo 
for (i in seq_along(marker_genes)) {
  
  gene <- marker_genes[i]

  expression_non_low <- marker_group_df[[gene]][
    marker_group_df$claudin_low_gallo == FALSE]
  expression_claudin_low <- marker_group_df[[gene]][
    marker_group_df$claudin_low_gallo == TRUE]
  
  marker_group_results$median_non_low[i] <- median(expression_non_low)
  marker_group_results$IQR_non_low[i] <- IQR(expression_non_low)
  
  marker_group_results$median_claudin_low[i] <- median(expression_claudin_low)
  marker_group_results$IQR_claudin_low[i] <- IQR(expression_claudin_low)
  
  marker_group_results$median_difference[i] <- median(expression_claudin_low) - 
    median(expression_non_low)
  
  comparison_test <- wilcox.test(x = expression_claudin_low, 
                                 y = expression_non_low,
                                 alternative = "two.sided", exact = FALSE)
  
  marker_group_results$p_value[i] <- comparison_test$p.value
}

marker_group_results$FDR <- p.adjust(marker_group_results$p_value, method = "BH")
marker_group_results

# si median difference negativo indica una menor expresión en claudin-low.
# si median difference positivo indica una mayor expresión en claudin-low.

# ******************************************************************************
# Visual of results
# ******************************************************************************
# Gallo_score distribution between claudin-low vs no-low groups
score_plot_df <- validation_df
score_plot_df$group <- factor(score_plot_df$claudin_low_gallo,
                              levels = c(FALSE, TRUE),
                              labels = c("No claudin-low", "Claudin-low"))
group_n <- table(score_plot_df$group)
# etiqueta grupos claudin
group_labels <- c(
  "No claudin-low" = paste0("No claudin-low\n(n = ", group_n["No claudin-low"], ")"),
  "Claudin-low" = paste0("Claudin-low\n(n = ", group_n["Claudin-low"], ")"))

wilcox_gallo <- wilcox.test(gallo_score ~ claudin_low_gallo, 
                            data = validation_df, exact = FALSE)
# etiqueta test wilcoxon
p_wilcox_label <- paste0("Wilcoxon: p ", format.pval(wilcox_gallo$p.value, digits = 2))

# anotación vertical
y_annotation <- max(score_plot_df$gallo_score, na.rm = TRUE) + 0.25
plot_gallo_score <- ggplot(score_plot_df, aes(x = group, y = gallo_score, fill = group)) +
  geom_violin(width = 0.75, trim = FALSE, alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.30) +
  annotate("text", x = 1.5, y = y_annotation, label = p_wilcox_label, size = 3.8) +
  scale_x_discrete(labels = group_labels) + 
  scale_fill_manual(values = c("No claudin-low" = "#A7A9AC", "Claudin-low" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Gallo score (mediana de z-scores)", 
       title = "Distribución del Gallo score") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

plot_gallo_score

# Marker genes expression in claudin-low vs non-low groups
marker_plot_df <- data.frame(gene = factor(
  rep(marker_genes, each = nrow(marker_group_df)), levels = marker_genes),
  expression = unlist(marker_group_df[marker_genes], use.names = FALSE),
  group = rep(marker_group_df$claudin_low_gallo, times = length(marker_genes)))

marker_plot_df$group <- factor(marker_plot_df$group, levels = c(FALSE, TRUE),
                               labels = c("No claudin-low", "Claudin-low"))

fdr_text <- ifelse(marker_group_results$FDR == 0, "FDR < 2.2e-16",
                   paste0("FDR = ", formatC(marker_group_results$FDR, format = "e", digits = 1)))

facet_labels <- setNames(paste0(marker_group_results$gene, "\n", fdr_text), 
                         marker_group_results$gene)

plot_markers <- ggplot(marker_plot_df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(width = 0.60, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.65, alpha = 0.20) +
  facet_wrap(~ gene, nrow = 1, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = c("No claudin-low" = "#A7A9AC", "Claudin-low" = "#D55E00")) +
  labs(x = NULL, y = expression(log[2](TPM + 1)), 
       title = "Expresión de marcadores del fenotipo claudin-low") +
  theme_classic(base_size = 11) + 
  theme(legend.position = "none", plot.title = element_text(face = "bold"), 
        strip.background = element_blank(), strip.text = element_text(face = "bold"), 
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_markers

# ******************************************************************************
# how much does the claudin-low signal track sample composition?
# ******************************************************************************
# calcular stromal, immune scores y pureza mediante TIDYESTIMATE
# para ello es necesario un df con los datos de expresion en log2(TPM+1), los id
# de muestras y o bien el entrez id o el gene symbol. Tidyestimate no acepta enseml.

stopifnot(is.matrix(tpm), is.numeric(tpm), all(c("gene_id", "gene_symbol") %in% 
                                                 colnames(tpm_df)))
# posiciones
gene_position <- match(rownames(tpm), tpm_df$gene_id)
if (anyNA(gene_position)) {
  missing_gene_ids <- rownames(tpm)[is.na(gene_position)]
  stop("No se han encontrado en tpm_df ", length(missing_gene_ids),
       " gene_id. Registros: ", paste(head(missing_gene_ids, 10), collapse = ", "))}
# gene map
gene_map <- tpm_df[gene_position, c("gene_id", "gene_symbol"), drop = FALSE]
gene_map$gene_id <- as.character(gene_map$gene_id)
gene_map$gene_symbol <- trimws(as.character(gene_map$gene_symbol))

duplicate_counts <- table(gene_map$gene_symbol)
if ((sum(duplicate_counts > 1)) >= 1) {
  stop("Hay ", (sum(duplicate_counts > 1)), " símbolos duplicados almenos 1 vez")
}

# corrección duplicados
valid_symbol <- (!is.na(gene_map$gene_symbol) & nzchar(gene_map$gene_symbol))
# se cogen los datos de tpm por simbolos que cumplen las condiciones normales y
# se suman las rows que tengan símbolos iguales manteniendo 1 símbolo
tpm_symbol <- rowsum(tpm[valid_symbol, , drop = FALSE],
  group = gene_map$gene_symbol[valid_symbol])

stopifnot(!anyDuplicated(rownames(tpm_symbol)), identical(colnames(tpm_symbol), colnames(tpm)))
stopifnot(nrow(tpm_symbol) == length(unique(gene_map$gene_symbol)))
tpm_log <- log2(tpm_symbol + 1)

# obtención de scores
tpm_estimate <- tpm_log |>
  filter_common_genes(id = "hgnc_symbol", tidy = FALSE, tell_missing = TRUE, 
                      find_alias = FALSE)
scores <- tpm_estimate |>
  estimate_score(is_affymetrix = FALSE)

# inserción en annotation
score_position <- match(annotation$rna_aliquot_id, scores$sample)
stopifnot(!anyNA(score_position))

annotation$stromal_score <- scores$stromal[score_position]
annotation$immune_score <- scores$immune[score_position]
annotation$estimate_score <- scores$estimate[score_position]

# cor Spearman gallo_score vs stromal/immune/estimate scores
score_names <- c("stromal_score", "immune_score", "estimate_score")
scores_long <- annotation |>
  select(gallo_score, all_of(score_names)) |> # selecciona las 4 columnas
  pivot_longer(cols = all_of(score_names), # pivot longer las pone una bajo la otra
               # y genera 2 columnas asociadas a los nombres y valores
               names_to = "score",
               values_to = "value") |>
  filter(is.finite(gallo_score), is.finite(value)) |>
  mutate(score = factor(score, levels = score_names))

score_correlations <- lapply(score_names, function(score_name) {
  df <- filter(scores_long, score == score_name) # se hace un df del panel
  # se calcula cor spearman de la pareja del panel
  test <- cor.test(df$gallo_score, df$value, method = "spearman", exact = FALSE)
  data.frame(score = score_name, n = nrow(df), rho = unname(test$estimate), 
             p_value = test$p.value)}) |>
  bind_rows() |>
  mutate(FDR = p.adjust(p_value, method = "BH"))
score_correlations

# se hacen las etiquetas con un df nuevo basado en score_correlations al que se
# le añade labels
plot_labels <- score_correlations |>
  mutate(label = paste0("n = ", n, "\nrho = ", round(rho, 3), 
                        "\nFDR ", format.pval(FDR, digits = 2)))

ggplot(scores_long, aes(x = gallo_score, y = value)) + # variables de los ejes
  geom_point(alpha = 0.4, size = 1.3) + # 1 punto por muestra
  # se añade una curva de tendencia
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#2166AC") +
  # se crean los paneles por separado en 1 fila
  facet_wrap(~ score, scales = "free_y", nrow = 1) +
  geom_text(data = plot_labels, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3.5) +
  labs(x = "Gallo score", y = "ESTIMATE-derived score") +
  theme_classic(base_size = 12)

# 2. Comparación claudin-low frente a non-low 
# Compararía cada score mediante Wilcoxon–Mann–Whitney bilateral para muestras independientes.



# 3. Ajuste e interpretación
# Aplicaría Benjamini–Hochberg por separado a las dos familias de preguntas: 
# los tres p-valores de correlación y los tres de comparación entre grupos. 
# Así quedarían seis pruebas, organizadas en dos bloques predefinidos.

# ...?


# testers
stopifnot(identical(colnames(counts), names(clusters)))
if (any(!is.finite(gene_sd) | gene_sd == 0)) {
  stop("Hay genes de la firma sin variabilidad entre muestras.")}



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