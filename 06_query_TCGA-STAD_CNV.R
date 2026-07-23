# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query_cnv <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  experimental.strategy = "Genotyping Array",
  platform = "Affymetrix SNP 6.0",
  workflow.type = "ASCAT3",
  access = "open"
)

GDCdownload(
  query_cnv,
  directory = tcga_dir
)

se_cnv <- GDCprepare(query_cnv,
                     directory = tcga_dir,
                     summarizedExperiment = TRUE)

saveRDS(
  se_cnv,
  file = file.path(tcga_dir, "Prepared", "TCGA_STAD_cnv_se.rds")
)
