# Install packages -------------------------------------------------------------
bioc_packages <- c(
  "BiocManager",
  "cBioPortalData"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# libraries --------------------------------------------------------------------
library(cBioPortalData)

# directory setup --------------------------------------------------------------
cbio_dir <- Sys.getenv("CBIO_DATA")

# query building ---------------------------------------------------------------
studies <- getStudies(cBioPortal()) # retrieve all studies from CBIOPORTAL

# study and id search
studies[
  grep("stomach|gastric", studies$name, ignore.case = TRUE),
  c("studyId", "name")
]

study_ids <- c("stad_tcga_pub")

# query download ---------------------------------------------------------------
for (id in study_ids) {
  downloadStudy(
    cancer_study_id = id,
    use_cache = cbio_dir,
    force = FALSE,
    ask = FALSE
  )
}
