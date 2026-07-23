# Install packages -------------------------------------------------------------
bioc_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# libraries --------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# directory setup --------------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# query building ---------------------------------------------------------------
query_mRNA <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  access = "open")

# query download ---------------------------------------------------------------
GDCdownload(
  query_mRNA,
  directory = tcga_dir)

se_mRNA <- GDCprepare(query_mRNA, 
                 directory = tcga_dir, 
                 summarizedExperiment = TRUE)

saveRDS(se_mRNA,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds")
)
