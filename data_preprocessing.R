####################
### Introduction ###
####################
# This script is used for performing the data pre-processing of in vivo mouse scRNA-seq data sets
# from two studies:
#   - De Micheli et al. (2020): GSE143435 and GSE143437 data sets
#   - Dell'Orso et al. (2019): GSE126834 data set

# The data sets are available in NCBI GEO:
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126834
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143435
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143437

# Data sets:
#  - GSE126834: 10X - total muscle wt, homeostatic MuSCs, Inj. 60h MuSCs, primary MB
#  - GSE143435: 10X - FACS sorted samples (d0, d2, d5, d7)
#  - GSE143437: 10X - non-FACS sorted samples

#####################
### Load packages ###
#####################
library(scater)
library(scran)
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)

# Set working dir to the source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################
### Utility functions ###
#########################

# Function to manipulate data so that the first column is moved to the rownames
firstColumnToRowNames <- function(data) {
  rownames(data) <- data[, 1]
  data <- data[, -1]
  return(data)
}

#####################
### Load data ###
#####################
# De Micheli et al. (2020): GSE143435 and GSE143437 data sets
# GSE143437 (non-FACS sorted samples: d0, d2, d5, d7)
GSE143437_meta <- read.delim("./../data/GSE143437/GSE143437_DeMicheli_MuSCatlas_metadata.txt")
GSE143437_raw <- read.delim("./../data/GSE143437/GSE143437_DeMicheli_MuSCatlas_rawdata.txt")
GSE143437_meta <- firstColumnToRowNames(GSE143437_meta)
GSE143437_raw <- firstColumnToRowNames(GSE143437_raw)

# GSE143435 (FACS sorted samples: d0, d2, d5, d7)
GSE143435_meta_d0 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D0_FACSatlas_metadata.txt")
GSE143435_raw_d0 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D0_FACSatlas_rawdata.txt")
GSE143435_meta_d2 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D2_FACSatlas_metadata.txt")
GSE143435_raw_d2 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D2_FACSatlas_rawdata.txt")
GSE143435_meta_d5 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D5_FACSatlas_metadata.txt")
GSE143435_raw_d5 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D5_FACSatlas_rawdata.txt")
GSE143435_meta_d7 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D7_FACSatlas_metadata.txt")
GSE143435_raw_d7 <- read.delim("./../data/GSE143435/GSE143435_DeMicheli_D7_FACSatlas_rawdata.txt")

GSE143435_raw_d0 <- firstColumnToRowNames(GSE143435_raw_d0)
GSE143435_meta_d0 <- firstColumnToRowNames(GSE143435_meta_d0)
GSE143435_raw_d2 <- firstColumnToRowNames(GSE143435_raw_d2)
GSE143435_meta_d2 <- firstColumnToRowNames(GSE143435_meta_d2)
GSE143435_raw_d5 <- firstColumnToRowNames(GSE143435_raw_d5)
GSE143435_meta_d5 <- firstColumnToRowNames(GSE143435_meta_d5)
GSE143435_raw_d7 <- firstColumnToRowNames(GSE143435_raw_d7)
GSE143435_meta_d7 <- firstColumnToRowNames(GSE143435_meta_d7)

library(dplyr) # Remove columns that don't exist in metadata
GSE143435_raw_d0 <- GSE143435_raw_d0 %>% select(rownames(GSE143435_meta_d0))
GSE143435_raw_d2 <- GSE143435_raw_d2 %>% select(rownames(GSE143435_meta_d2))
GSE143435_raw_d5 <- GSE143435_raw_d5 %>% select(rownames(GSE143435_meta_d5))
GSE143435_raw_d7 <- GSE143435_raw_d7 %>% select(rownames(GSE143435_meta_d7))

# Dell'Orso et al. (2019): GSE126834 data set
samples <- c("./../data/GSE126834/homeostatic_muscs_1", 
             "./../data/GSE126834/inj_60h_muscs_1", 
             "./../data/GSE126834/total_muscle_wt_1",
             "./../data/GSE126834/homeostatic_muscs_2",
             "./../data/GSE126834/inj_60h_muscs_2",
             "./../data/GSE126834/primary_MB",
             "./../data/GSE126834/total_muscle_wt_2")

GSE126834 <- Read10X(data.dir = samples)

#########################
### Build SCE objects ###
#########################
# GSE143437 (non-FACS)
sce_GSE143437 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143437_raw)),
                                      colData = GSE143437_meta)

# GSE143435 (FACS: d0, d2, d5, d7)
sce_GSE143435_d0 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d0)),
                                         colData = GSE143435_meta_d0)
sce_GSE143435_d2 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d2)),
                                         colData = GSE143435_meta_d2)
sce_GSE143435_d5 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d5)),
                                         colData = GSE143435_meta_d5)
sce_GSE143435_d7 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d7)),
                                         colData = GSE143435_meta_d7)
dim(GSE143435_raw_d0)
dim(GSE143435_meta_d0)

# GSE126834
sce_GSE126834 <- SingleCellExperiment(assays = list(counts = GSE126834))

# Save objects in RDS to avoid long computational times
saveRDS(sce_GSE143437, file = "sce_GSE143437")
saveRDS(sce_GSE143435_d0, file = "sce_GSE143435_d0")
saveRDS(sce_GSE143435_d2, file = "sce_GSE143435_d2")
saveRDS(sce_GSE143435_d5, file = "sce_GSE143435_d5")
saveRDS(sce_GSE143435_d7, file = "sce_GSE143435_d7")
saveRDS(sce_GSE126834, file = "sce_GSE126834")