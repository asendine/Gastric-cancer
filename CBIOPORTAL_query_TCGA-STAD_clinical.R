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

# STUDY SEARCH -----------------------------------------------------------------
studies <- getStudies(cBioPortal()) # retrieve all studies from CBIOPORTAL

# study and id search
studies[
  grep("stomach", studies$name, ignore.case = TRUE),
  c("studyId", "name")
]

#-------------------------------------------------------------------------------
# STUDY 1 DATA -----------------------------------------------------------------
#-------------------------------------------------------------------------------

study_id <- c("stad_tcga_pub") # se determina el dataset n = ~250

# Retreiving cBioPortal data ---------------------------------------------------
# idealmente se haría un multiassay:
stad_tcga_pub_ma <- cBioDataPack(
  cancer_study_id = "stad_tcga_pub",
  use_cache = cbio_dir,
  cleanup = TRUE,
  ask = FALSE
)

# pero si no funciona porque hay variables que no casan, se sigue descargando todo
# manualmente

# query download ---------------------------------------------------------------
downloadStudy( 
  cancer_study_id = study_id,
  use_cache = cbio_dir,
  force = FALSE,
  ask = FALSE
  )

# Al descomprimir el archivo se obtienen diferentes datos. Para saber qué datos 
# es de interés importar se obtienen los nombres de todas las columnas:

# primero se descomprime manualmente el tar.gz y se define stad_dir
stad_dir <- file.path(cbio_dir, "stad_tcga_pub", "stad_tcga_pub")

# CHECK COLUMNAS ---------------------------------------------------------------
# all file routes are retrieved from the decompressed original file 
files <- list.files(
  stad_dir,
  pattern = "\\.(txt|seg)$",
  full.names = TRUE
)

# se obtiene un objeto con los archivos y sus columnas
columns <- lapply(files, function(file) {
  names(read.delim(
    file,
    comment.char = "#",
    nrows = 0,
    check.names = FALSE
  ))
})

names(columns) <- basename(files)
columns

# en principio se puede ver que se va a necesitar: data_clinical_patient, 
# data_clinical_sample, data_cna (copy number alterations), data_linear_cna, 
# data_mrna_seq_v2_rsem, data_mutations.

# OBTENCIÓN DATOS --------------------------------------------------------------
clin_pat_pub <- read.delim(file.path(stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)

clin_sample_pub <- read.delim(file.path(stad_dir, "data_clinical_sample.txt"),
                              comment.char = "#",check.names = FALSE)

cna_pub <- read.delim(file.path(stad_dir, "data_cna.txt"),check.names = FALSE)

cna_lin_pub <- read.delim(file.path(stad_dir, "data_linear_cna.txt"),
                         check.names = FALSE)

mrna_pub <- read.delim(file.path(stad_dir, "data_mrna_seq_v2_rsem.txt"),
                   check.names = FALSE)

mut_pub <- read.delim(file.path(stad_dir, "data_mutations.txt"),comment.char = "#",
                  check.names = FALSE)


#-------------------------------------------------------------------------------
# STUDY 2 DATA -----------------------------------------------------------------
#-------------------------------------------------------------------------------

study_id <- c("stad_tcga_gdc") # se determina el dataset n = ~450

# Retreiving cBioPortal data ---------------------------------------------------
# idealmente se haría un multiassay:
stad_tcga_gdc_ma <- cBioDataPack(
  cancer_study_id = "stad_tcga_gdc",
  use_cache = cbio_dir,
  cleanup = TRUE,
  ask = FALSE
)

# pero si no funciona porque hay variables que no casan, se sigue descargando todo
# manualmente

# query download ---------------------------------------------------------------
downloadStudy( 
  cancer_study_id = study_id,
  use_cache = cbio_dir,
  force = FALSE,
  ask = FALSE
)

# Al descomprimir el archivo se obtienen diferentes datos. Para saber qué datos 
# es de interés importar se obtienen los nombres de todas las columnas:

# primero se descomprime manualmente el tar.gz y se define stad_dir
stad_dir <- file.path(cbio_dir, "stad_tcga_gdc", "stad_tcga_gdc")

# CHECK COLUMNAS ---------------------------------------------------------------
# all file routes are retrieved from the decompressed original file 
files <- list.files(
  stad_dir,
  pattern = "\\.(txt|seg)$",
  full.names = TRUE
)

# se obtiene un objeto con los archivos y sus columnas
columns <- lapply(files, function(file) {
  names(read.delim(
    file,
    comment.char = "#",
    nrows = 0,
    check.names = FALSE
  ))
})

names(columns) <- basename(files)
columns

# en principio se puede ver que se va a necesitar: data_clinical_patient, 
# data_clinical_sample, data_cna (copy number alterations), data_linear_cna, 
# data_mrna_seq_v2_rsem, data_mutations.

# OBTENCIÓN DATOS --------------------------------------------------------------
clin_pat_gdc <- read.delim(file.path(stad_dir, "data_clinical_patient.txt"),
                           comment.char = "#",check.names = FALSE)

clin_sample_gdc <- read.delim(file.path(stad_dir, "data_clinical_sample.txt"),
                              comment.char = "#",check.names = FALSE)

cna_gdc <- read.delim(file.path(stad_dir, "data_cna.txt"), check.names = FALSE)

cna_lin_gdc <- read.delim(file.path(stad_dir, "data_linear_cna.txt"),
                          check.names = FALSE)

mrna_gdc <- read.delim(file.path(stad_dir, "data_mrna_seq_v2_rsem.txt"),
                       check.names = FALSE)

mut_gdc <- read.delim(file.path(stad_dir, "data_mutations.txt"), comment.char = "#",
                      check.names = FALSE)