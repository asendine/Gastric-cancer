# Install packages CBIO --------------------------------------------------------
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

# directory definition ---------------------------------------------------------
tcga_dir <- Sys.getenv("TCGA_DATA")

# building a query -------------------------------------------------------------
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  access = "open")

GDCdownload(
  query,
  directory = tcga_dir)

se <- GDCprepare(query)

saveRDS(se,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_rnaseq_se.rds")
)