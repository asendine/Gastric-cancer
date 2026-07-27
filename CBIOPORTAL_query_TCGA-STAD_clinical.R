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
# si interesa obtener varios estudios se puede hacer un for pero no es típico
downloadStudy( 
  cancer_study_id = id,
  use_cache = cbio_dir,
  force = FALSE,
  ask = FALSE
  )

# Retreiving cBioPortal data ---------------------------------------------------
# idealmente se haría un multiassay:

#stad_tcga_pub <- cBioDataPack(
#  cancer_study_id = "stad_tcga_pub",
#  use_cache = cbio_dir,
#  cleanup = TRUE,
#  ask = FALSE
#)

# pero no funciona, da error con el fichero de stad_tcga_pub.
# Queda la opción de ir obteniendo los archivos uno a uno. Al descomprimir el
# archivo se obtienen diferentes datos. Para saber qué datos es de interés importar
# se obtienen los nombres de todas las columnas:

# se redefine stad_dir ya que se ha generado una carpeta extra al descomprimir:
stad_dir <- file.path(cbio_dir, "stad_tcga_pub", "stad_tcga_pub")

# CHECK COLUMNAS ---------------------------------------------------------------
# all file rutes are retrieved from the decompresssed original file 
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
clinical_patient <- read.delim(file.path(stad_dir, "data_clinical_patient.txt"),
                               comment.char = "#",check.names = FALSE)

clinical_sample <- read.delim(file.path(stad_dir, "data_clinical_sample.txt"),
                              comment.char = "#",check.names = FALSE)

cna <- read.delim(file.path(stad_dir, "data_cna.txt"),check.names = FALSE)

cna_linear <- read.delim(file.path(stad_dir, "data_linear_cna.txt"),
                         check.names = FALSE)

mrna <- read.delim(file.path(stad_dir, "data_mrna_seq_v2_rsem.txt"),
                   check.names = FALSE)

mut <- read.delim(file.path(stad_dir, "data_mutations.txt"),comment.char = "#",
                  check.names = FALSE)

