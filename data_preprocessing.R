####################
### Introduction ###
####################
# This script is used for performing the data pre-processing of in vivo mouse scRNA-seq data sets from two studies:
#   - De Micheli et al. (2020): GSE143435 and GSE143437 data sets
#   - Dell'Orso et al. (2019): GSE126834 data set

# The data sets are available in NCBI GEO:
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126834
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143435
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143437

# Set working dir to the source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#####################
### Load packages ###
#####################
library(scater)
library(scran)
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)
library(dplyr)
library(docstring)
library(ggplot2)

#########################
### Utility functions ###
#########################

firstColumnToRowNames <- function(data) {
  #' Function to manipulate data so that the first column is moved to the rownames
  rownames(data) <- data[, 1]
  data <- data[, -1]
  return(data)
}


removeGeneExpLessThanX <- function(data, x) {
  #' Remove genes (from a sce object) that are expressed in less than x cells
  keep_feature <- rowSums(counts(data) > 0) >= x
  data <- data[keep_feature, ]
  return(data)
}


performCellQC <- function(sce) {
  #' Add per cell QC metrics and detect mitochondrial transcripts
  location <- rowRanges(sce)
  is.mito <- any(seqnames(location) == "MT")
  sce <- addPerCellQC(sce, subsets = list(Mito = is.mito))
  return(sce)
}


findOutliers <- function(sce) {
  #' Detect outliers based on median absolute deviation (MAD)
  # Library size outliers
  qc.lib.sce <- isOutlier(sce$sum, log = TRUE, type = "both", nmads = 4)
  # Expressed genes outliers
  qc.nexprs.sce <- isOutlier(sce$detected, log = TRUE, type = "both", nmads = 4)
  # Mitochondrial content outliers
  qc.mito.sce <- sce$subsets_Mito_percent > 20 # or isOutlier(sce$subsets_Mito_percent, type = "higher")
  # All outliers
  discard <- qc.lib.sce | qc.nexprs.sce | qc.mito.sce
  # Summarize the number of cells removed for each reason
  print(data.frame(LibSize = sum(qc.lib.sce),
                   NExprs = sum(qc.nexprs.sce),
                   MitoProp = sum(qc.mito.sce),
                   Total = sum(discard)))
  
  return(discard)
}


checkQCMetricCorrelation <- function(sce) {
  #' Check if qc metrics correlate
  cdf <- colData(sce)[c("subsets_Mito_percent", "sum", "detected")]
  cors <- cor(as.matrix(cdf))
  # Print the name of the input object and the correlation matrix
  print(substitute(sce))
  return(cors)
}


performNormalization <- function(sce) {
  #' Computes library size factors, sum factors, and deconvolution size factors and log-transforms the data
  # Compute library size factors
  lib.sf <- librarySizeFactors(sce)
  # Perform clustering
  set.seed(100)
  quickie <- quickCluster(sce)
  print(table(quickie))
  # compute sum factors
  sce <- computeSumFactors(sce, cluster = quickie)
  # Compute deconvolution size factors
  deconv.sf <- sizeFactors(sce)
  # plot lib sf vs. deconv sf.
  plot(x = lib.sf,
       y = deconv.sf,
       xlab = "library size factor",
       ylab = "deconvolution size factor",
       log = "xy",
       pch = 16)
  abline(a=0, b=1, col="red")
  # perform log normalization
  sce <- logNormCounts(sce)
}


checkDiagnosticPlots <- function(sce) {
  #' Create a diagnostic plot of total counts vs. mito%
  return(
  gridExtra::grid.arrange(
    plotColData(sce,
                x = "sum",
                y = "subsets_Mito_percent",
                colour_by = "discard")
    )
  )
}

checkForEmptyDroplets <- function(sce) {
  #' Checks for empty droplets and subsets sce to retain only detected cells
  # 'NOTE: # CellRanger v3 automatically performs cell calling so there is no need for calling emptyDrops()
  
  set.seed(100)
  # Distinguish between droplets containing cells and ambient RNA
  e.out <- emptyDrops(counts(sce),
                      lower = 100,
                      test.ambient = TRUE)
  # Summarize the results
  print(substitute(sce))
  print(summary(e.out))
  # Subset sce object to retain only the detected cells
  sce <- sce[, which(e.out$FDR <= 0.001)]
  return(sce)
}


#####################
### Load data ###
#####################

# --------------------------------------------------------------- #
# De Micheli et al. (2020) 
# --------------------------------------------------------------- #
# Data set consists of non-FACS (GSE143437) and FACS (GSE143435) sorted samples.
# Features: 
# - Injury: [Day 0, Day 2, Day 5, Day7]
# - SampleID: FACS: [D0_FACS, D2_FACS, D5_FACS, D7_FACS]
#             non-FACS: [D0_A, D0_B, D0_Cv3, D2_C, D2_D, D5_A, D5_B,
#                        D5_C, D7_C, D7_D]
# - cell_annotation: ["Anti-flammatory machophages", "FAPs", "B cells",
#                     "MuSCs and progenitors", "B cells", "NK cells", ...]
# --------------------------------------------------------------- #

# Load GSE143437 data set (non-FACS sorted samples: d0, d2, d5, d7)
GSE143437_meta <- read.delim("./../data/GSE143437/GSE143437_DeMicheli_MuSCatlas_metadata.txt")
GSE143437_raw <- read.delim("./../data/GSE143437/GSE143437_DeMicheli_MuSCatlas_rawdata.txt")
GSE143437_meta <- firstColumnToRowNames(GSE143437_meta)
GSE143437_raw <- firstColumnToRowNames(GSE143437_raw)

# Load GSE143435 data set (FACS sorted samples: d0, d2, d5, d7)
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

# Remove columns that do not exist in metadata
GSE143435_raw_d0 <- GSE143435_raw_d0 %>% select(rownames(GSE143435_meta_d0))
GSE143435_raw_d2 <- GSE143435_raw_d2 %>% select(rownames(GSE143435_meta_d2))
GSE143435_raw_d5 <- GSE143435_raw_d5 %>% select(rownames(GSE143435_meta_d5))
GSE143435_raw_d7 <- GSE143435_raw_d7 %>% select(rownames(GSE143435_meta_d7))


# --------------------------------------------------------------- #
# Dell'Orso et al. (2019)
# --------------------------------------------------------------- #
# Data set consisting of 7 samples
# - "homeostatic_MuSCs_rep1 / rep2"
# - "inj_60h_MuSCs_rep1 / rep 2"
# - "Primary_MB"
# - "total muscle rep1 / rep2"
# --------------------------------------------------------------- #

# Get the sample directory paths
samples <- c("./../data/GSE126834/homeostatic_muscs_1", 
             "./../data/GSE126834/homeostatic_muscs_2",
             "./../data/GSE126834/inj_60h_muscs_1", 
             "./../data/GSE126834/inj_60h_muscs_2",
             "./../data/GSE126834/primary_MB",
             "./../data/GSE126834/total_muscle_wt_1",
             "./../data/GSE126834/total_muscle_wt_2")

# Read samples separately
GSE126834_hom1 <- Read10X(data.dir = samples[1])
GSE126834_hom2 <- Read10X(data.dir = samples[2])
GSE126834_inj1 <- Read10X(data.dir = samples[3])
GSE126834_inj2 <- Read10X(data.dir = samples[4])
GSE126834_mb <- Read10X(data.dir = samples[5])
GSE126834_tot1 <- Read10X(data.dir = samples[6])
GSE126834_tot2 <- Read10X(data.dir = samples[6])

#########################
### Build SCE objects ###
#########################
# --------------------------------------------------------------- #
# De Micheli et al. (2020): GSE143435 and GSE143437 data sets
# --------------------------------------------------------------- #
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

# Combine data
sce_GSE143435_d0$FACS <- TRUE
sce_GSE143435_d2$FACS <- TRUE
sce_GSE143435_d5$FACS <- TRUE
sce_GSE143435_d7$FACS <- TRUE
sce_GSE143437$FACS <- FALSE

sce.deMicheli <- list(FACS_d0 = sce_GSE143435_d0,
                      FACS_d2 = sce_GSE143435_d2,
                      FACS_d5 = sce_GSE143435_d5,
                      FACS_d7 = sce_GSE143435_d7,
                      non_FACS = sce_GSE143437)

# Save sce object
saveRDS(sce.deMicheli, file = "sce_deMicheli")

# free memory
remove(sce_GSE143435_d0)
remove(sce_GSE143435_d2)
remove(sce_GSE143435_d5)
remove(sce_GSE143435_d7)
remove(sce_GSE143437)
remove(GSE143435_raw_d0)
remove(GSE143435_raw_d2)
remove(GSE143435_raw_d5)
remove(GSE143435_raw_d7)
remove(GSE143437_raw)

# --------------------------------------------------------------- #
# Dell'Orso et al. (2019): GSE126834 data set
# --------------------------------------------------------------- #
# Create sce objects
sce_GSE126834_hom1 <- SingleCellExperiment(assays = list(counts = GSE126834_hom1))
sce_GSE126834_hom2 <- SingleCellExperiment(assays = list(counts = GSE126834_hom2))
sce_GSE126834_inj1 <- SingleCellExperiment(assays = list(counts = GSE126834_inj1))
sce_GSE126834_inj2 <- SingleCellExperiment(assays = list(counts = GSE126834_inj2))
sce_GSE126834_mb <- SingleCellExperiment(assays = list(counts = GSE126834_mb))
sce_GSE126834_tot1 <- SingleCellExperiment(assays = list(counts = GSE126834_tot1))
sce_GSE126834_tot2 <- SingleCellExperiment(assays = list(counts = GSE126834_tot2))


# Assign sample groups to the sce objects
sce_GSE126834_hom1$sample <- "homeostatic_MuSCs_rep1"
sce_GSE126834_hom2$sample <- "homeostatic_MuSCs_rep2"
sce_GSE126834_inj1$sample <- "inj_60h_MuSCs_rep1"
sce_GSE126834_inj2$sample <- "inj_60h_MuSCs_rep2"
sce_GSE126834_mb$sample <- "Primary_MB"
sce_GSE126834_tot1$sample <- "total_muscle_wt_rep1"
sce_GSE126834_tot2$sample <- "total_muscle_wt_rep2"

# Merge the samples into a single sce object (samples separated by the $sample column)
sce.dellOrso <- cbind(sce_GSE126834_hom1, sce_GSE126834_hom2,
                      sce_GSE126834_inj1, sce_GSE126834_inj2,
                      sce_GSE126834_mb, sce_GSE126834_tot1,
                      sce_GSE126834_tot2)

# Save sce object
saveRDS(sce.dellOrso, file = "sce_dellOrso")

# free memory
remove(sce_GSE126834_hom1)
remove(sce_GSE126834_hom2)
remove(sce_GSE126834_inj1)
remove(sce_GSE126834_inj2)
remove(sce_GSE126834_mb)
remove(sce_GSE126834_tot1)
remove(sce_GSE126834_tot2)
remove(GSE126834_hom1)
remove(GSE126834_hom2)
remove(GSE126834_inj1)
remove(GSE126834_inj2)
remove(GSE126834_mb)
remove(GSE126834_tot1)
remove(GSE126834_tot2)

#######################
### Quality control ###
#######################
# --------------------------------------------------------------- #
# Notes about quality control
# --------------------------------------------------------------- #
# Motivation: remove low-quality libraries (cells) from the data

# Metrics of low quality:
#  - Cells with small library sizes (Library size = total sum of counts across all relevant features for each cell)
#  - Cells with only few expressed genes (number of genes with non-zero counts for that cell)
#  - Proportion of reads mapped to mitochondrial genome (high proportions = poor quality)
#  - (Proportion of reads mapped to spike-in transcripts (high proportions = poor quality))

# In practice:
#  - perCellQualityMetrics() calculates all the above metrics
#     * sum = total count for each cell
#     * detected = # of detected genes
#     * subsets_mito_percent = % of reads mapped to mitochondrial transcripts
#     * altexps_ERCC_percent = % of reads mapped to ERCC transcripts (spike-ins)
#  - addPerCellQC() does the same but just appends te per-cell QC metrics to the colData of sce object

# Identifying outliers
#  - Outliers can be detected based on median absolute deviation (MAD) from the median value
#  - A value is considered an outlier if it is more than 3 MADs from the median in the "problematic" direction
# --------------------------------------------------------------- #

# Remove genes that are expressed in less than 3 cells
sce.dellOrso <- removeGeneExpLessThanX(sce.dellOrso, 3)
for (n in names(sce.deMicheli)) {
  sce.deMicheli[[n]] <- removeGeneExpLessThanX(sce.deMicheli[[n]], 3)
  }

# Compute per cell quality control metrics
sce.dellOrso <- performCellQC(sce.dellOrso)
for (n in names(sce.deMicheli)) {
  sce.deMicheli[[n]] <- performCellQC(sce.deMicheli[[n]])
  }

# Save objects
saveRDS(sce.dellOrso, file = "sce_dellOrso_qc")
saveRDS(sce.deMicheli, file = "sce_dellMicheli_qc")

# Identify outliers (low quality cells)
sce.dellOrso$discard <- findOutliers(sce.dellOrso)

for (n in names(sce.deMicheli)) {
  sce.deMicheli[[n]]$discard <- findOutliers(sce.deMicheli[[n]])
  }

# Summarize outliers
summary(sce.dellOrso$discard)

for (n in names(sce.deMicheli)) {
  print(n)
  print(summary(sce.deMicheli[[n]]$discard))
}

# Discard outliers
sce.dellOrso.f <- sce.dellOrso[, sce.dellOrso$discard == FALSE]

sce.deMicheli.f <- sce.deMicheli
for (n in names(sce.deMicheli.f)) {
  sce.deMicheli.f[[n]] <- sce.deMicheli.f[[n]][, sce.deMicheli.f[[n]]$discard == FALSE]
  }

# Save objects
saveRDS(sce.deMicheli.f, file = "sce_deMicheli_f")
saveRDS(sce.dellOrso.f, file = "sce_dellOrso_f")


#####################
### Normalization ###
#####################
# --------------------------------------------------------------- #
# Notes about normalization
# --------------------------------------------------------------- #
# Motivation: There are systematic differences in seq. coverage between libraries. These differences arise from technical differences.
#             Normalization aims to remove these differences such that they do not interfere with comparisons of the expression profiles between cells

# In practice:
#  - Scaling normalization: dividing all counts for each cell by a cell-specific scaling factor, often called a “size factor”
#     * The size factor for each cell represents the estimate of the relative bias in that cell, so division of its counts by its size factor should remove that bias.
#     * librarySizeFactors() function
#     * Once size factors are computed, we use logNormCounts() to compute normalized exp. values for each cell
#         + Divides each gene with the appropriate size factor
#         + Log-transforms the normalized values (creating a new assay "logcounts")
#         + -> log-values represent log-fold changes in expression
#  - Normalization by deconvolution: normalize on summed expression values from pools of cells
#     * quickCluster() + calculateSumFactors()
# --------------------------------------------------------------- #

performNormalization(sce.dellOrso.f)
for(n in names(sce.deMicheli.f)) {
  performNormalization(sce.deMicheli.f[[n]])
}

# Cell cycle scoring (https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html)












