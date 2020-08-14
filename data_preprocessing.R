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

# Remove genes from a sce object that are expressed in less than x cells
removeGeneExpLessThanX <- function(data, x) {
  keep_feature <- rowSums(counts(data) > 0) >= x
  data <- data[keep_feature, ]
  return(data)
}

# Perform pel cell QC and detect mitochondrial transcripts
mitoPerCellQC <- function(sce) {
  location <- rowRanges(sce)
  is.mito <- any(seqnames(location) == "MT")
  sce <- addPerCellQC(sce, subsets = list(Mito = is.mito))
  return(sce)
}

# Detect outliers (sum, detected, mito_percent)
findOutliers <- function(sce) {
  qc.lib.sce <- isOutlier(sce$sum, log = TRUE, type = "both", nmads = 4)
  qc.nexprs.sce <- isOutlier(sce$detected, log = TRUE, type = "both", nmads = 4)
  qc.mito.sce <- sce$subsets_Mito_percent > 20
  
  discard <- qc.lib.sce | qc.nexprs.sce | qc.mito.sce
  return(discard)
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

# Save objects to avoid long computational times
saveRDS(sce_GSE143437, file = "sce_GSE143437")
saveRDS(sce_GSE143435_d0, file = "sce_GSE143435_d0")
saveRDS(sce_GSE143435_d2, file = "sce_GSE143435_d2")
saveRDS(sce_GSE143435_d5, file = "sce_GSE143435_d5")
saveRDS(sce_GSE143435_d7, file = "sce_GSE143435_d7")
saveRDS(sce_GSE126834, file = "sce_GSE126834")

#######################
### Quality control ###
#######################
# Load files from a file (optional)
#sce_GSE143437 <- readRDS("sce_GSE143437")
#sce_GSE143435_d0 <- readRDS("sce_GSE143435_d0")
#sce_GSE143435_d2 <- readRDS("sce_GSE143435_d2")
#sce_GSE143435_d5 <- readRDS("sce_GSE143435_d5")
#sce_GSE143435_d7 <- readRDS("sce_GSE143435_d7")
#sce_GSE126834 <- readRDS("sce_GSE126834")

# Remove genes that are expressed in less than 3 cells
sce_GSE143437 <- removeGeneExpLessThanX(sce_GSE143437, 3)
sce_GSE143435_d0 <- removeGeneExpLessThanX(sce_GSE143435_d0, 3)
sce_GSE143435_d2 <- removeGeneExpLessThanX(sce_GSE143435_d2, 3)
sce_GSE143435_d5 <- removeGeneExpLessThanX(sce_GSE143435_d5, 3)
sce_GSE143435_d7 <- removeGeneExpLessThanX(sce_GSE143435_d7, 3)
sce_GSE126834 <- removeGeneExpLessThanX(sce_GSE126834, 3)

# Retrieve mitochondrial transcripts and perform per cell QC
sce_GSE143437 <- mitoPerCellQC(sce_GSE143437)
sce_GSE143435_d0 <- mitoPerCellQC(sce_GSE143435_d0)
sce_GSE143435_d2 <- mitoPerCellQC(sce_GSE143435_d2)
sce_GSE143435_d5 <- mitoPerCellQC(sce_GSE143435_d5)
sce_GSE143435_d7 <- mitoPerCellQC(sce_GSE143435_d7)
sce_GSE126834 <- mitoPerCellQC(sce_GSE126834)

# Save objects
saveRDS(sce_GSE143437, file = "sce_GSE143437_qc")
saveRDS(sce_GSE143435_d0, file = "sce_GSE143435_d0_qc")
saveRDS(sce_GSE143435_d2, file = "sce_GSE143435_d2_qc")
saveRDS(sce_GSE143435_d5, file = "sce_GSE143435_d5_qc")
saveRDS(sce_GSE143435_d7, file = "sce_GSE143435_d7_qc")
saveRDS(sce_GSE126834, file = "sce_GSE126834_qc")

# Identify outliers (low quality cells)
discard_sce_GSE143437 <- findOutliers(sce_GSE143437)
discard_sce_GSE143435_d0 <- findOutliers(sce_GSE143435_d0)
discard_sce_GSE143435_d2 <- findOutliers(sce_GSE143435_d2)
discard_sce_GSE143435_d5 <- findOutliers(sce_GSE143435_d5)
discard_sce_GSE143435_d7 <- findOutliers(sce_GSE143435_d7)
discard_sce_GSE126834 <- findOutliers(sce_GSE126834)

sce_GSE143437$discard <- discard_sce_GSE143437
sce_GSE143435_d0$discard <- discard_sce_GSE143435_d0
sce_GSE143435_d2$discard <- discard_sce_GSE143435_d2
sce_GSE143435_d5$discard <- discard_sce_GSE143435_d5
sce_GSE143435_d7$discard <- discard_sce_GSE143435_d7
sce_GSE126834$discard <- discard_sce_GSE126834

# Discard outliers
sce_GSE143437.f <- sce_GSE143437[, sce_GSE143437$discard == FALSE]
sce_GSE143435_d0.f <- sce_GSE143435_d0[, sce_GSE143435_d0$discard == FALSE]
sce_GSE143435_d2.f <- sce_GSE143435_d2[, sce_GSE143435_d2$discard == FALSE]
sce_GSE143435_d5.f <- sce_GSE143435_d5[, sce_GSE143435_d5$discard == FALSE]
sce_GSE143435_d7.f <- sce_GSE143435_d7[, sce_GSE143435_d7$discard == FALSE]
sce_GSE126834.f <- sce_GSE126834[, sce_GSE126834$discard == FALSE]

# Save objects
saveRDS(sce_GSE143437.f, file = "sce_GSE143437_qc_f")
saveRDS(sce_GSE143435_d0.f, file = "sce_GSE143435_d0_qc_f")
saveRDS(sce_GSE143435_d2.f, file = "sce_GSE143435_d2_qc_f")
saveRDS(sce_GSE143435_d5.f, file = "sce_GSE143435_d5_qc_f")
saveRDS(sce_GSE143435_d7.f, file = "sce_GSE143435_d7_qc_f")
saveRDS(sce_GSE126834.f, file = "sce_GSE126834_qc_f")

#sce_GSE143437.f <- readRDS("sce_GSE143437_qc_f")
#sce_GSE143435_d0.f <- readRDS("sce_GSE143435_d0_qc_f")
#sce_GSE143435_d2.f <- readRDS("sce_GSE143435_d2_qc_f")
#sce_GSE143435_d5.f <- readRDS("sce_GSE143435_d5_qc_f")
#sce_GSE143435_d7.f <- readRDS("sce_GSE143435_d7_qc_f")
#sce_GSE126834.f <- readRDS("sce_GSE126834_qc_f")

#####################
### Normalization ###
#####################
lib.sf <- librarySizeFactors(sce_GSE143437.f)
set.seed(100)
quickie <- quickCluster(sce_GSE143437.f)
table(quickie)



# Log normalization


# Cell cycle scoring (https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html)












