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
library(dplyr)
library(docstring)

# Set working dir to the source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
  #' Remove genes from a sce object that are expressed in less than x cells
  keep_feature <- rowSums(counts(data) > 0) >= x
  data <- data[keep_feature, ]
  return(data)
}


performQC <- function(sce) {
  #' Perform pel cell QC and detect mitochondrial transcripts
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
  data.frame(LibSize = sum(qc.lib.sce),
             NExprs = sum(qc.nexprs.sce),
             MitoProp = sum(qc.mito.sce),
             Total = sum(discard))
  
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

# Read samples separately (homeostatic MUSCs, injured 60h MUSCs, Primary Myoblasts)
# TODO: Find out if total muscle should also be added to the data set
GSE126834_hom1 <- Read10X(data.dir = samples[1])
GSE126834_hom2 <- Read10X(data.dir = samples[2])
GSE126834_inj1 <- Read10X(data.dir = samples[3])
GSE126834_inj2 <- Read10X(data.dir = samples[4])
GSE126834_mb <- Read10X(data.dir = samples[5])


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

# TODO: Find out how and if the datasets should be merged together prior to normalization?
#sce_de_micheli <- cbind(sce_GSE143435_d0, sce_GSE143435_d2, sce_GSE143435_d5, sce_GSE143435_d7, sce_GSE143437)

# --------------------------------------------------------------- #
# Dell'Orso et al. (2019): GSE126834 data set
# --------------------------------------------------------------- #
# Create sce objects
sce_GSE126834_hom1 <- SingleCellExperiment(assays = list(counts = GSE126834_hom1))
sce_GSE126834_hom2 <- SingleCellExperiment(assays = list(counts = GSE126834_hom2))
sce_GSE126834_inj1 <- SingleCellExperiment(assays = list(counts = GSE126834_inj1))
sce_GSE126834_inj2 <- SingleCellExperiment(assays = list(counts = GSE126834_inj2))
sce_GSE126834_mb <- SingleCellExperiment(assays = list(counts = GSE126834_mb))

# Assign sample groups to the sce objects
sce_GSE126834_hom1$sample <- "homeostatic_MuSCs_rep1"
sce_GSE126834_hom2$sample <- "homeostatic_MuSCs_rep2"
sce_GSE126834_inj1$sample <- "inj_60h_MuSCs_rep1"
sce_GSE126834_inj2$sample <- "inj_60h_MuSCs_rep2"
sce_GSE126834_mb$sample <- "Primary_MB"

# Merge the samples into a single sce object (samples separated by the $sample column)
sce_GSE126834 <- cbind(sce_GSE126834_hom1, sce_GSE126834_hom2,
                       sce_GSE126834_inj1, sce_GSE126834_inj2,
                       sce_GSE126834_mb)



# Save sce objects
saveRDS(sce_GSE143437, file = "sce_GSE143437")
saveRDS(sce_GSE143435_d0, file = "sce_GSE143435_d0")
saveRDS(sce_GSE143435_d2, file = "sce_GSE143435_d2")
saveRDS(sce_GSE143435_d5, file = "sce_GSE143435_d5")
saveRDS(sce_GSE143435_d7, file = "sce_GSE143435_d7")
saveRDS(sce_GSE126834_merged, file = "sce_GSE126834_merged")



#######################
### Quality control ###
#######################

#### Notes on QC ######
## Motivation: remove low-quality libraries (cells) from the data ##

## Metrics of low quality: ##
#   - Cells with small library sizes (Library size = total sum of counts across all relevant features for each cell)
#   - Cells with only few expressed genes (number of genes with non-zero counts for that cell)
#   - Proportion of reads mapped to mitochondrial genome (high proportions = poor quality)
#   - (Proportion of reads mapped to spike-in transcripts (high proportions = poor quality))

## In practice:
#   - perCellQualityMetrics() calculates all the above metrics
#       * sum = total count for each cell
#       * detected = # of detected genes
#       * subsets_mito_percent = % of reads mapped to mitochondrial transcripts
#       * altexps_ERCC_percent = % of reads mapped to ERCC transcripts (spike-ins)
#  - addPerCellQC() does the same but just appends te per-cell QC metrics to the colData of sce object

# Identifying outliers
#  - Outliers can be detected based on median absolute deviation (MAD) from the median value
#  - A value is considered an outlier if it is more than 3 MADs from the median in the "problematic" direction
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
sce_GSE143437 <- performQC(sce_GSE143437)
sce_GSE143435_d0 <- performQC(sce_GSE143435_d0)
sce_GSE143435_d2 <- performQC(sce_GSE143435_d2)
sce_GSE143435_d5 <- performQC(sce_GSE143435_d5)
sce_GSE143435_d7 <- performQC(sce_GSE143435_d7)
sce_GSE126834 <- performQC(sce_GSE126834)

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

sce_GSE143437.f <- readRDS("sce_GSE143437_qc_f")
sce_GSE143435_d0.f <- readRDS("sce_GSE143435_d0_qc_f")
sce_GSE143435_d2.f <- readRDS("sce_GSE143435_d2_qc_f")
sce_GSE143435_d5.f <- readRDS("sce_GSE143435_d5_qc_f")
sce_GSE143435_d7.f <- readRDS("sce_GSE143435_d7_qc_f")
sce_GSE126834.f <- readRDS("sce_GSE126834_qc_f")

#######################
# QC diagnostic plots #
# Motivation: It is good practice to inspect the distributions of QC metrics to identify possible problems
#  - in ideal case we would see normal distributions that would justify the 3 MAD treshold used in outlier detection

# Check if QC metrics correlate in the data sets
checkQCMetricCorrelation(sce_GSE143437.f)
checkQCMetricCorrelation(sce_GSE143435_d0.f)
checkQCMetricCorrelation(sce_GSE143435_d2.f)
checkQCMetricCorrelation(sce_GSE143435_d5.f)
checkQCMetricCorrelation(sce_GSE143435_d7.f)
checkQCMetricCorrelation(sce_GSE126834.f)

# Check the relationship between total counts vs. mito%
checkDiagnosticPlots(sce_GSE143437.f)
checkDiagnosticPlots(sce_GSE143435_d0.f)
checkDiagnosticPlots(sce_GSE143435_d2.f)
checkDiagnosticPlots(sce_GSE143435_d5.f)
checkDiagnosticPlots(sce_GSE143435_d7.f)
checkDiagnosticPlots(sce_GSE126834.f)


#######################

#######################
# Cell calling
# - Motivation: Call cells from empty droplets based on the observed expression profiles.
# - In practice: emptyDrops() function tests whether the expression profile for each cell barcode is significantly different from the ambient RNA pool
#    * This allows us to discriminate between well-sequenced empty droplets and droplets derived from cells with little RNA

# !! Move to the utility functions

checkForEmptyDroplets <- function(sce) {
  # NOTE: # CellRanger v3 automatically performs cell calling so there is no need for calling emptyDrops()
  set.seed(100)
  # Distinguish between droplets containing cells and ambient RNA
  e.out <- emptyDrops(counts(sce),
                      lower = 100,
                      test.ambient = TRUE)
  # Summarize the results
  summary(e.out)
  # Subset sce object to retain only the detected cells
  sce <- sce[, which(e.out$FDR <= 0.001)]
  return(sce)
}
#######################

#####################
### Normalization ###
#####################
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

performNormalization(sce_GSE143437.f)
performNormalization(sce_GSE143435_d0.f)
performNormalization(sce_GSE143435_d2.f)
performNormalization(sce_GSE143435_d5.f)
performNormalization(sce_GSE143435_d7.f)
performNormalization(sce_GSE126834.f)

# Cell cycle scoring (https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html)












