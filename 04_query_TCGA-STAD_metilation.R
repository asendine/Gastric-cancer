# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  access = "open"
)

GDCdownload(
  query,
  directory = tcga_dir)

se <- GDCprepare(query)

saveRDS(se,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_methylation_450K_se.rds")
)