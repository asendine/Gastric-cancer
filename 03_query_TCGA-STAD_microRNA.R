# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query_microRNA <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(
  query_microRNA,
  directory = tcga_dir)

se_microRNA <- GDCprepare(query_microRNA, 
                          directory = tcga_dir, 
                          summarizedExperiment = TRUE)

saveRDS(se_microRNA,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_mirna_se.rds")
)
