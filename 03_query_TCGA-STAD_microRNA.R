# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(
  query,
  directory = tcga_dir)

se <- GDCprepare(query)

saveRDS(se,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_mirna_se.rds")
)