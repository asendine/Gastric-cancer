# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number Scores"
)

GDCdownload(
  query,
  directory = tcga_dir
)

cnv_gene <- GDCprepare(query)

saveRDS(
  cnv_gene,
  file = file.path(tcga_dir, "Prepared", "TCGA_STAD_gene_cnv.rds")
)