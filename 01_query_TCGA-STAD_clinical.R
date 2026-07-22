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
results <- GDCquery_clinic(project = "TCGA-STAD", 
                           type = "clinical")

# saving the data --------------------------------------------------------------
saveRDS(results,
        file = file.path(tcga_dir, "Prepared", "TCGA_STAD_clinical.rds")
)
