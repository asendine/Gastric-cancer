# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

GDCdownload(
  query,
  directory = tcga_dir
)

maf <- GDCprepare(query)

saveRDS(
  maf,
  file = file.path(tcga_dir, "Prepared", "TCGA_STAD_maf.rds")
)