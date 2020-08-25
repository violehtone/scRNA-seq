####################
### Introduction ###
####################
# This script is used for performing the data pre-processing of 3 scRNA-seq data sets.
# One dataset is primary data and the other two are used as reference datasets.
# Reference data sets:
#   - De Micheli et al. (2020): GSE143435 data set
#   - Dell'Orso et al. (2019): GSE126834 data set

# The data sets are available in NCBI GEO:
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126834
#   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143435

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
library(edgeR)

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

#####################
### Load data ###
#####################

# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #
counts <- Read10X("./../data/filtered_feature_bc_matrix/")
counts <- counts[rowSums(counts) > 0, ]

cellcodes <- as.data.frame(counts@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes
cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes))
samples<-c("control","cap50","cap50_r4h","cap50_r8h","cap50_r16h")
cellcodes$samples <- as.vector(samples[cellcodes$libcodes])

# --------------------------------------------------------------- #
# De Micheli et al. (2020) 
# --------------------------------------------------------------- #
# Data set consisting of 4 samples: d0, d2, d5, d7
# Features:
# - Injury: [Day 0, Day 2, Day 5, Day7]
# - sampleID: [D0_FACS, D2_FACS, D5_FACS, D7_FACS]
# - cell_annotation: ["Anti-flammatory machophages", "FAPs", "B cells",
#                     "MuSCs and progenitors", "B cells", "NK cells", ...]
# --------------------------------------------------------------- #

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

# Remove cells that do not exist in metadata
GSE143435_raw_d0 <- GSE143435_raw_d0 %>% select(rownames(GSE143435_meta_d0))
GSE143435_raw_d2 <- GSE143435_raw_d2 %>% select(rownames(GSE143435_meta_d2))
GSE143435_raw_d5 <- GSE143435_raw_d5 %>% select(rownames(GSE143435_meta_d5))
GSE143435_raw_d7 <- GSE143435_raw_d7 %>% select(rownames(GSE143435_meta_d7))

# Add zero counts to make the dimensions (genes) of data sets similar
# d0
GSE143435_raw_d0[setdiff(rownames(GSE143435_raw_d2), rownames(GSE143435_raw_d0)), ] <- 0
GSE143435_raw_d0[setdiff(rownames(GSE143435_raw_d5), rownames(GSE143435_raw_d0)), ] <- 0
GSE143435_raw_d0[setdiff(rownames(GSE143435_raw_d7), rownames(GSE143435_raw_d0)), ] <- 0
# d2
GSE143435_raw_d2[setdiff(rownames(GSE143435_raw_d0), rownames(GSE143435_raw_d2)), ] <- 0
GSE143435_raw_d2[setdiff(rownames(GSE143435_raw_d5), rownames(GSE143435_raw_d2)), ] <- 0
GSE143435_raw_d2[setdiff(rownames(GSE143435_raw_d7), rownames(GSE143435_raw_d2)), ] <- 0
# d5
GSE143435_raw_d5[setdiff(rownames(GSE143435_raw_d0), rownames(GSE143435_raw_d5)), ] <- 0
GSE143435_raw_d5[setdiff(rownames(GSE143435_raw_d2), rownames(GSE143435_raw_d5)), ] <- 0
GSE143435_raw_d5[setdiff(rownames(GSE143435_raw_d7), rownames(GSE143435_raw_d5)), ] <- 0
# d7
GSE143435_raw_d7[setdiff(rownames(GSE143435_raw_d0), rownames(GSE143435_raw_d7)), ] <- 0
GSE143435_raw_d7[setdiff(rownames(GSE143435_raw_d2), rownames(GSE143435_raw_d7)), ] <- 0
GSE143435_raw_d7[setdiff(rownames(GSE143435_raw_d5), rownames(GSE143435_raw_d7)), ] <- 0

# Fix colnames of d5 metadata (for some reason it had percent.mito instead of percent_mito)
names(GSE143435_meta_d5)[names(GSE143435_meta_d5) == "percent.mito"] <- "percent_mito"

# Add unique identifier for column names of each dataset
colnames(GSE143435_raw_d0) <- paste(colnames(GSE143435_raw_d0), "d0", sep = "_")
colnames(GSE143435_raw_d2) <- paste(colnames(GSE143435_raw_d2), "d2", sep = "_")
colnames(GSE143435_raw_d5) <- paste(colnames(GSE143435_raw_d5), "d5", sep = "_")
colnames(GSE143435_raw_d7) <- paste(colnames(GSE143435_raw_d7), "d7", sep = "_")

rownames(GSE143435_meta_d0) <- paste(rownames(GSE143435_meta_d0), "d0", sep = "_")
rownames(GSE143435_meta_d2) <- paste(rownames(GSE143435_meta_d2), "d2", sep = "_")
rownames(GSE143435_meta_d5) <- paste(rownames(GSE143435_meta_d5), "d5", sep = "_")
rownames(GSE143435_meta_d7) <- paste(rownames(GSE143435_meta_d7), "d7", sep = "_")

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
             "./../data/GSE126834/primary_MB")

# Read samples separately
GSE126834_hom1 <- Read10X(data.dir = samples[1])
GSE126834_hom2 <- Read10X(data.dir = samples[2])
GSE126834_inj1 <- Read10X(data.dir = samples[3])
GSE126834_inj2 <- Read10X(data.dir = samples[4])
GSE126834_mb <- Read10X(data.dir = samples[5])

# Add unique identifier for column names of each dataset
colnames(GSE126834_hom1) <- paste(colnames(GSE126834_hom1), "hom1", sep = "_")
colnames(GSE126834_hom2) <- paste(colnames(GSE126834_hom2), "hom2", sep = "_")
colnames(GSE126834_inj1) <- paste(colnames(GSE126834_inj1), "inj1", sep = "_")
colnames(GSE126834_inj2) <- paste(colnames(GSE126834_inj2), "inj2", sep = "_")
colnames(GSE126834_mb) <- paste(colnames(GSE126834_mb), "mb", sep = "_")


#########################
### Build SCE objects ###
#########################
# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #
sce <- SingleCellExperiment(assays=list(counts=counts))
colData(sce)$sample <- as.vector(cellcodes$samples)

# Save sce object
saveRDS(sce, file = "sce")

# Free space
remove(counts)

# --------------------------------------------------------------- #
# De Micheli et al. (2020): GSE143435 data set
# --------------------------------------------------------------- #
# merge raw data into a single object
sce_GSE143435_d0 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d0)),
                                         colData = GSE143435_meta_d0)
sce_GSE143435_d2 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d2)),
                                         colData = GSE143435_meta_d2)
sce_GSE143435_d5 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d5)),
                                         colData = GSE143435_meta_d5)
sce_GSE143435_d7 <- SingleCellExperiment(assays = list(counts = as.matrix(GSE143435_raw_d7)),
                                         colData = GSE143435_meta_d7)

# Merge samples into a single sce object
sce.deMicheli <- cbind(sce_GSE143435_d0, sce_GSE143435_d2,
                       sce_GSE143435_d5, sce_GSE143435_d7,
                       deparse.level = 1)

# Save sce object
saveRDS(sce.deMicheli, file = "sce_deMicheli")

# free memory
remove(sce_GSE143435_d0)
remove(sce_GSE143435_d2)
remove(sce_GSE143435_d5)
remove(sce_GSE143435_d7)
remove(GSE143435_raw_d0)
remove(GSE143435_raw_d2)
remove(GSE143435_raw_d5)
remove(GSE143435_raw_d7)
remove(GSE143435_meta_d0)
remove(GSE143435_meta_d2)
remove(GSE143435_meta_d5)
remove(GSE143435_meta_d7)

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
sce.dellOrso <- cbind(sce_GSE126834_hom1, sce_GSE126834_hom2,
                      sce_GSE126834_inj1, sce_GSE126834_inj2,
                      sce_GSE126834_mb,
                      deparse.level = 1)

# Save sce object
saveRDS(sce.dellOrso, file = "sce_dellOrso")

# free memory
remove(sce_GSE126834_hom1)
remove(sce_GSE126834_hom2)
remove(sce_GSE126834_inj1)
remove(sce_GSE126834_inj2)
remove(sce_GSE126834_mb)
remove(GSE126834_hom1)
remove(GSE126834_hom2)
remove(GSE126834_inj1)
remove(GSE126834_inj2)
remove(GSE126834_mb)

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
#  - addPerCellQC() does the same but just appends the QC metrics to the colData of sce object

# Identifying outliers
#  - Outliers can be detected based on median absolute deviation (MAD) from the median value
#  - A value is considered an outlier if it is more than 3 MADs from the median in the "problematic" direction
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #
# Compute per cell quality control metrics
sce <- addPerCellQC(sce,subsets=list(Mito=grep("mt-", rownames(sce))))
rowData(sce)$mito <- FALSE
rowData(sce)$mito[grep("^mt-",rownames(sce))] <- TRUE

# Find outliers from total counts for each cell (sum)
qc.lib.control <- isOutlier(sce$sum[sce$sample=="control"], log=TRUE, type="both",nmads =4)
summary(qc.lib.control)

qc.lib.cap50 <- isOutlier(sce$sum[sce$sample=="cap50"], log=TRUE, type="both",nmads =4)
summary(qc.lib.cap50)

qc.lib.cap50_r4h <- isOutlier(sce$sum[sce$sample=="cap50_r4h"], log=TRUE, type="both",nmads =4)
summary(qc.lib.cap50_r4h)

qc.lib.cap50_r8h <- isOutlier(sce$sum[sce$sample=="cap50_r8h"], log=TRUE, type="both",nmads =4)
summary(qc.lib.cap50_r8h)

qc.lib.cap50_r16h <- isOutlier(sce$sum[sce$sample=="cap50_r16h"], log=TRUE, type="both",nmads =4)
summary(qc.lib.cap50_r16h)

# Find outliers from total number of detected genes
qc.nexprs.control <- isOutlier(sce$detected[sce$sample=="control"], log=TRUE, type="both",nmads =4)
summary(qc.nexprs.control)

qc.nexprs.cap50 <- isOutlier(sce$detected[sce$sample=="cap50"], log=TRUE, type="both",nmads =4)
summary(qc.nexprs.cap50)

qc.nexprs.cap50_r4h <- isOutlier(sce$detected[sce$sample=="cap50_r4h"], log=TRUE, type="both",nmads =4)
summary(qc.nexprs.cap50_r4h)

qc.nexprs.cap50_r8h <- isOutlier(sce$detected[sce$sample=="cap50_r8h"], log=TRUE, type="both",nmads =4)
summary(qc.nexprs.cap50_r8h)

qc.nexprs.cap50_r16h <- isOutlier(sce$detected[sce$sample=="cap50_r16h"], log=TRUE, type="both",nmads =4)
summary(qc.nexprs.cap50_r16h)

# Find outliers from percentage of reads mapped to mitochondrial transcripts
qc.mito.control <- sce$subsets_Mito_percent[sce$sample=="control"] > 20
summary(qc.mito.control)

qc.mito.cap50 <- sce$subsets_Mito_percent[sce$sample=="cap50"] > 20
summary(qc.mito.cap50)

qc.mito.cap50_r4h <- sce$subsets_Mito_percent[sce$sample=="cap50_r4h"] > 20
summary(qc.mito.cap50_r4h)

qc.mito.cap50_r8h <- sce$subsets_Mito_percent[sce$sample=="cap50_r8h"] > 20
summary(qc.mito.cap50_r8h)

qc.mito.cap50_r16h <- sce$subsets_Mito_percent[sce$sample=="cap50_r16h"] > 20

# Check if qc metrics correlate

cdf <- colData(sce)[c("subsets_Mito_percent","sum","detected")]
cors <- cor(as.matrix(cdf))
cors

# Combine discards

discard.control <- qc.lib.control | qc.nexprs.control | qc.mito.control
discard.cap50 <- qc.lib.cap50 | qc.nexprs.cap50 | qc.mito.cap50
discard.cap50_r4h <- qc.lib.cap50_r4h | qc.nexprs.cap50_r4h | qc.mito.cap50_r4h
discard.cap50_r8h <- qc.lib.cap50_r8h | qc.nexprs.cap50_r8h | qc.mito.cap50_r8h
discard.cap50_r16h <- qc.lib.cap50_r16h | qc.nexprs.cap50_r16h | qc.mito.cap50_r16h


# Summarize the number of cells removed for each reason
discardSummary.control <- DataFrame(LibSize=sum(qc.lib.control), NExprs=sum(qc.nexprs.control),
                                    MitoProp=sum(qc.mito.control), Total=sum(discard.control))

discardSummary.cap50 <- DataFrame(LibSize=sum(qc.lib.cap50), NExprs=sum(qc.nexprs.cap50),
                                  MitoProp=sum(qc.mito.cap50), Total=sum(discard.cap50))

discardSummary.cap50_r4h <- DataFrame(LibSize=sum(qc.lib.cap50_r4h), NExprs=sum(qc.nexprs.cap50_r4h),
                                      MitoProp=sum(qc.mito.cap50_r4h), Total=sum(discard.cap50_r4h))

discardSummary.cap50_r8h <- DataFrame(LibSize=sum(qc.lib.cap50_r8h), NExprs=sum(qc.nexprs.cap50_r8h),
                                      MitoProp=sum(qc.mito.cap50_r8h), Total=sum(discard.cap50_r8h))

discardSummary.cap50_r16h <- DataFrame(LibSize=sum(qc.lib.cap50_r16h), NExprs=sum(qc.nexprs.cap50_r16h),
                                       MitoProp=sum(qc.mito.cap50_r16h), Total=sum(discard.cap50_r16h))

discardSummary.control
discardSummary.cap50
discardSummary.cap50_r4h
discardSummary.cap50_r8h
discardSummary.cap50_r16h

sce$discard <- c(discard.control,discard.cap50,discard.cap50_r4h,discard.cap50_r8h,discard.cap50_r16h)

# Plot QC metrics
gridExtra::grid.arrange(
  plotColData(sce, x="sample", y="sum",colour_by = "discard",point_size=0.3) +
    scale_y_log10() + ggtitle("Total counts") + guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce, x="sample", y="detected",colour_by = "discard",point_size=0.3) + 
    scale_y_log10() + ggtitle("Detected genes") + guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce, x="sample", y="subsets_Mito_percent",colour_by = "discard",point_size=0.3) + 
    ggtitle("Percent mitochondrial reads") + guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce, x="sample", y="percent_top_500",colour_by = "discard",point_size=0.3) + 
    ggtitle("Percent top 500 genes") + guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  ncol=1
)

# Filter outliers

sce.f <- sce[,sce$discard==FALSE]

# Check for cell type enrichment in the discarded pool

lost <- calculateAverage(counts(sce)[,c(discard.control,discard.cap50,discard.cap50_r4h,discard.cap50_r8h,discard.cap50_r16h)])
kept <- calculateAverage(counts(sce)[,!c(discard.control,discard.cap50,discard.cap50_r4h,discard.cap50_r8h,discard.cap50_r16h)])

library(edgeR)

logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Plot gene logFC between discarded and kept cells (blue = mitochondrial genes)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[rowData(sce)$mito==TRUE], logFC[rowData(sce)$mito==TRUE], col="dodgerblue", pch=16)

changed <- logFC>0.5

# Up in discarded
rownames(sce)[changed]

# Save filtered object
saveRDS(sce.f,file="sce_f")

# --------------------------------------------------------------- #
# De Micheli et al. data set
# --------------------------------------------------------------- #
# sce.deMicheli <- readRDS("sce_deMicheli")

# Compute per cell quality control metrics
sce.deMicheli <- addPerCellQC(sce.deMicheli, subsets = list(Mito=grep("mt-", rownames(sce.deMicheli))))
rowData(sce.deMicheli)$mito <- FALSE
rowData(sce.deMicheli)$mito[grep("^mt-", rownames(sce.deMicheli))] <- TRUE

# Find outliers from total counts for each cell (sum)
qc.lib.d0 <- isOutlier(sce.deMicheli$sum[sce.deMicheli$sampleID == "D0_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.d0)
qc.lib.d2 <- isOutlier(sce.deMicheli$sum[sce.deMicheli$sampleID == "D2_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.d2)
qc.lib.d5 <- isOutlier(sce.deMicheli$sum[sce.deMicheli$sampleID == "D5_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.d5)
qc.lib.d7 <- isOutlier(sce.deMicheli$sum[sce.deMicheli$sampleID == "D7_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.d7)

# Find outliers from total number of detected genes
qc.nexprs.d0 <- isOutlier(sce.deMicheli$detected[sce.deMicheli$sampleID == "D0_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.nexprs.d0)
qc.nexprs.d2 <- isOutlier(sce.deMicheli$detected[sce.deMicheli$sampleID == "D2_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.nexprs.d2)
qc.nexprs.d5 <- isOutlier(sce.deMicheli$detected[sce.deMicheli$sampleID == "D5_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.nexprs.d5)
qc.nexprs.d7 <- isOutlier(sce.deMicheli$detected[sce.deMicheli$sampleID == "D7_FACS"], log = TRUE, type = "both", nmads = 3)
summary(qc.nexprs.d7)

# Find outliers from percentage of reads mapped to mitochondrial transcripts
qc.mito.d0 <- sce.deMicheli$subsets_Mito_percent[sce.deMicheli$sampleID == "D0_FACS"] > 20
summary(qc.mito.d0)
qc.mito.d2 <- sce.deMicheli$subsets_Mito_percent[sce.deMicheli$sampleID == "D2_FACS"] > 20
summary(qc.mito.d2)
qc.mito.d5 <- sce.deMicheli$subsets_Mito_percent[sce.deMicheli$sampleID == "D5_FACS"] > 20
summary(qc.mito.d5)
qc.mito.d7 <- sce.deMicheli$subsets_Mito_percent[sce.deMicheli$sampleID == "D7_FACS"] > 20
summary(qc.mito.d7)

# Check if QC metrics correlate
cdf.deMicheli <- colData(sce.deMicheli)[c("subsets_Mito_percent","sum","detected")]
cors.deMicheli <- cor(as.matrix(cdf.deMicheli))
cors.deMicheli

# Combine discards
discard.d0 <- qc.lib.d0 | qc.nexprs.d0 | qc.mito.d0
discard.d2 <- qc.lib.d2 | qc.nexprs.d2 | qc.mito.d2
discard.d5 <- qc.lib.d5 | qc.nexprs.d5 | qc.mito.d5
discard.d7 <- qc.lib.d7 | qc.nexprs.d7 | qc.mito.d7

# Summarize the number of cells removed for each reason
discardSummary.d0 <- DataFrame(LibSize = sum(qc.lib.d0),
                               NExprs = sum(qc.nexprs.d0),
                               MitoProp = sum(qc.mito.d0),
                               Total = sum(discard.d0))
discardSummary.d2 <- DataFrame(LibSize = sum(qc.lib.d2),
                               NExprs = sum(qc.nexprs.d2),
                               MitoProp = sum(qc.mito.d2),
                               Total = sum(discard.d2))
discardSummary.d5 <- DataFrame(LibSize = sum(qc.lib.d5),
                               NExprs = sum(qc.nexprs.d5),
                               MitoProp = sum(qc.mito.d5),
                               Total = sum(discard.d5))
discardSummary.d7 <- DataFrame(LibSize = sum(qc.lib.d7),
                               NExprs = sum(qc.nexprs.d7),
                               MitoProp = sum(qc.mito.d7),
                               Total = sum(discard.d7))

discardSummary.d0
discardSummary.d2
discardSummary.d5
discardSummary.d7

sce.deMicheli$discard <- c(discard.d0, discard.d2, discard.d5, discard.d7)

# Plot Quality control metrics
#TODO: find out why sce.deMicheli$subsets_Mito_percent[sce.deMicheli$sampleID == "D5_FACS"] returns only zeros
gridExtra::grid.arrange(
  plotColData(sce.deMicheli, x = "sampleID", y = "sum", colour_by = "discard", point_size = 0.3) +
    scale_y_log10() +
    ggtitle("Total counts") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.deMicheli, x = "sampleID", y = "detected", colour_by = "discard", point_size = 0.3) +
    scale_y_log10() +
    ggtitle("Total counts") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.deMicheli, x = "sampleID", y = "subsets_Mito_percent", colour_by = "discard", point_size = 0.3) +
    ggtitle("Total counts") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.deMicheli, x="sampleID", y="percent_top_500",colour_by = "discard",point_size=0.3) + 
    ggtitle("Percent top 500 genes") + guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  ncol = 1
)

# Filter outliers
sce.deMicheli.f <- sce.deMicheli[, sce.deMicheli$discard == FALSE]

# Check for cell type enrichments in the discarded pool
lost <- calculateAverage(counts(sce.deMicheli)[,c(discard.d0, discard.d2, discard.d5, discard.d7)])
kept <- calculateAverage(counts(sce.deMicheli)[,!c(discard.d0, discard.d2, discard.d5, discard.d7)])

logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Plot gene logFC between discarded and kept cells (blue = mitochondrial genes)

plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[rowData(sce.deMicheli)$mito==TRUE], logFC[rowData(sce.deMicheli)$mito==TRUE], col="dodgerblue", pch=16)

changed <- logFC>0.5

# Up in discarded
rownames(sce.deMicheli)[changed]

# Save filtered sce file
saveRDS(sce.deMicheli.f, file = "sce_deMicheli_f")


# --------------------------------------------------------------- #
# Dell'Orso et al. data set
# --------------------------------------------------------------- #
#sce.dellOrso <- readRDS("sce_dellOrso")

# Compute per cell quality control metrics
sce.dellOrso <- addPerCellQC(sce.dellOrso, subsets = list(Mito=grep("mt-", rownames(sce.dellOrso))))
rowData(sce.dellOrso)$mito <- FALSE
rowData(sce.dellOrso)$mito[grep("^mt-", rownames(sce.dellOrso))] <- TRUE

# Find outliers from total counts for each cell (sum)
qc.lib.hom1 <- isOutlier(sce.dellOrso$sum[sce.dellOrso$sample == "homeostatic_MuSCs_rep1"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.hom1)
qc.lib.hom2 <- isOutlier(sce.dellOrso$sum[sce.dellOrso$sample == "homeostatic_MuSCs_rep2"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.hom2)
qc.lib.inj1 <- isOutlier(sce.dellOrso$sum[sce.dellOrso$sample == "inj_60h_MuSCs_rep1"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.inj1)
qc.lib.inj2 <- isOutlier(sce.dellOrso$sum[sce.dellOrso$sample == "inj_60h_MuSCs_rep2"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.inj2)
qc.lib.pmb <- isOutlier(sce.dellOrso$sum[sce.dellOrso$sample == "Primary_MB"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.pmb)

#Find outliers from total number of detected genes
qc.nexprs.hom1 <- isOutlier(sce.dellOrso$detected[sce.dellOrso$sample == "homeostatic_MuSCs_rep1"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.hom1)
qc.nexprs.hom2 <- isOutlier(sce.dellOrso$detected[sce.dellOrso$sample == "homeostatic_MuSCs_rep2"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.hom2)
qc.nexprs.inj1 <- isOutlier(sce.dellOrso$detected[sce.dellOrso$sample == "inj_60h_MuSCs_rep1"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.inj1)
qc.nexprs.inj2 <- isOutlier(sce.dellOrso$detected[sce.dellOrso$sample == "inj_60h_MuSCs_rep2"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.inj2)
qc.nexprs.pmb <- isOutlier(sce.dellOrso$detected[sce.dellOrso$sample == "Primary_MB"], log = TRUE, type = "both", nmads = 3)
summary(qc.lib.pmb)

# Find outliers from percentage of reads mapped to mitochondrial transcripts
qc.mito.hom1 <- sce.dellOrso$subsets_Mito_percent[sce.dellOrso$sample == "homeostatic_MuSCs_rep1"] > 20
summary(qc.mito.hom1)
qc.mito.hom2 <- sce.dellOrso$subsets_Mito_percent[sce.dellOrso$sample == "homeostatic_MuSCs_rep2"] > 20
summary(qc.mito.hom2)
qc.mito.inj1 <- sce.dellOrso$subsets_Mito_percent[sce.dellOrso$sample == "inj_60h_MuSCs_rep1"] > 20
summary(qc.mito.inj1)
qc.mito.inj2 <- sce.dellOrso$subsets_Mito_percent[sce.dellOrso$sample == "inj_60h_MuSCs_rep2"] > 20
summary(qc.mito.inj2)
qc.mito.pmb <- sce.dellOrso$subsets_Mito_percent[sce.dellOrso$sample == "Primary_MB"] > 20
summary(qc.lib.pmb)

# Check if QC metrics correlate
cdf.dellOrso <- colData(sce.dellOrso)[c("subsets_Mito_percent","sum","detected")]
cors.dellOrso <- cor(as.matrix(cdf.dellOrso))
cors.dellOrso

# Combine discards
discard.hom1 <- qc.lib.hom1 | qc.nexprs.hom1 | qc.mito.hom1
discard.hom2 <- qc.lib.hom2 | qc.nexprs.hom2 | qc.mito.hom2
discard.inj1 <- qc.lib.inj1 | qc.nexprs.inj1 | qc.mito.inj1
discard.inj2 <- qc.lib.inj2 | qc.nexprs.inj2 | qc.mito.inj2
discard.pmb <- qc.lib.pmb | qc.nexprs.pmb | qc.mito.pmb

# Summarize the number of cells removed for each reason
discardSummary.hom1 <- DataFrame(LibSize = sum(qc.lib.hom1),
                               NExprs = sum(qc.nexprs.hom1),
                               MitoProp = sum(qc.mito.hom1),
                               Total = sum(discard.hom1))
discardSummary.hom2 <- DataFrame(LibSize = sum(qc.lib.hom2),
                                 NExprs = sum(qc.nexprs.hom2),
                                 MitoProp = sum(qc.mito.hom2),
                                 Total = sum(discard.hom2))
discardSummary.inj1 <- DataFrame(LibSize = sum(qc.lib.inj1),
                                 NExprs = sum(qc.nexprs.inj1),
                                 MitoProp = sum(qc.mito.inj1),
                                 Total = sum(discard.inj1))
discardSummary.inj2 <- DataFrame(LibSize = sum(qc.lib.inj2),
                                 NExprs = sum(qc.nexprs.inj2),
                                 MitoProp = sum(qc.mito.inj2),
                                 Total = sum(discard.inj2))
discardSummary.pmb <- DataFrame(LibSize = sum(qc.lib.pmb),
                                 NExprs = sum(qc.nexprs.pmb),
                                 MitoProp = sum(qc.mito.pmb),
                                 Total = sum(discard.pmb))

discardSummary.hom1
discardSummary.hom2
discardSummary.inj1
discardSummary.inj2
discardSummary.pmb

sce.dellOrso$discard <- c(discard.hom1, discard.hom2, discard.inj1, discard.inj2, discard.pmb)

# Plot QC metrics
gridExtra::grid.arrange(
  plotColData(sce.dellOrso, x="sample", y="sum",colour_by = "discard",point_size=0.3) +
    scale_y_log10() + 
    ggtitle("Total counts") + 
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.dellOrso, x="sample", y="detected",colour_by = "discard",point_size=0.3) + 
    scale_y_log10() +
    ggtitle("Detected genes") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.dellOrso, x="sample", y="subsets_Mito_percent",colour_by = "discard",point_size=0.3) + 
    ggtitle("Percent mitochondrial reads") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  plotColData(sce.dellOrso, x="sample", y="percent_top_500",colour_by = "discard",point_size=0.3) + 
    ggtitle("Percent top 500 genes") +
    guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))),
  ncol=1
)


# Filter outliers
sce.dellOrso.f <- sce.dellOrso[, sce.dellOrso$discard == FALSE]

# Check for cell type enrichment in the discarded pool

lost <- calculateAverage(counts(sce.dellOrso)[,c(discard.hom1, discard.hom2, discard.inj1, discard.inj2, discard.pmb)])
kept <- calculateAverage(counts(sce.dellOrso)[,!c(discard.hom1, discard.hom2, discard.inj1, discard.inj2, discard.pmb)])

logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Plot gene logFC between discarded and kept cells (blue = mitochondrial genes)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[rowData(sce.dellOrso)$mito==TRUE], logFC[rowData(sce.dellOrso)$mito==TRUE], col="dodgerblue", pch=16)

changed <- logFC>0.5

# Up in discarded
rownames(sce.dellOrso)[changed]

# Save filtered sce file
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

# Make new sce objects for normalization purposes
sce.n <- sce.f
sce.dellOrso.n <- sce.dellOrso.f
sce.deMicheli.n <- sce.deMicheli.f

# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #
# Compute size factors
lib.sf <- librarySizeFactors(sce.n)
summary(lib.sf)

set.seed(100)
quickie <- quickCluster(sce.n,block=sce.n$sample)
table(quickie)

sce.n <- computeSumFactors(sce.n, cluster=quickie)

deconv.sf <- sizeFactors(sce.n)

# Plot library size factor vs. deconvolution size factor
plot(lib.sf, deconv.sf, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

# Compute log-normalized expression values for each cell
sce.n <- logNormCounts(sce.n)

# Save normalized data
saveRDS(sce.n, file = "sce_n")

# --------------------------------------------------------------- #
# Dell'Orso et al.
# --------------------------------------------------------------- #
# Compute size factors
lib.sf.dellOrso <- librarySizeFactors(sce.dellOrso.n)
summary(lib.sf.dellOrso)

set.seed(100)
clust.dellOrso <- quickCluster(sce.dellOrso.n)

sce.dellOrso.n <- computeSumFactors(sce.dellOrso.n,
                                    cluster = clust.dellOrso,
                                    min.mean = 0.1)

deconv.sf.dellOrso <- sizeFactors(sce.dellOrso.n)
summary(deconv.sf.dellOrso)

# Plot library size factor vs. deconvolution size factor
plot(lib.sf.dellOrso, deconv.sf.dellOrso, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

# Compute log-normalized expression values for each cell
sce.dellOrso.n <- logNormCounts(sce.dellOrso.n)

# Save normalized data
saveRDS(sce.dellOrso.n, file = "sce_dellOrso_n")

# --------------------------------------------------------------- #
# De Micheli et al.
# --------------------------------------------------------------- #
# Compute size factors
lib.sf.deMicheli <- librarySizeFactors(sce.deMicheli.n)
summary(lib.sf.deMicheli)

set.seed(100)
clust.deMicheli <- quickCluster(sce.deMicheli.n)

sce.deMicheli.n <- computeSumFactors(sce.deMicheli.n,
                                     cluster = clust.deMicheli,
                                     min.mean = 0.1)

deconv.sf.deMicheli <- sizeFactors(sce.deMicheli.n)
summary(deconv.sf.deMicheli)

# Plot library size factors vs. deconvolution size factor
plot(lib.sf.deMicheli, deconv.sf.deMicheli, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

# Compute log-normalized expression values for each cell
sce.deMicheli.n <- logNormCounts(sce.deMicheli.n)

# Save normalized data
saveRDS(sce.deMicheli.n, file = "sce_deMicheli_n")


#########################################
### Cell cycle scoring and regression ###
#########################################
# --------------------------------------------------------------- #
# Notes about cell cycle scoring and regression
# --------------------------------------------------------------- #
# Solution 1: 
# we assign each cell a score, based on its expression of G2/M and S phase markers. 
# Cells expressing neither are likely in G1 phase

# In practice:
#  - We assign scores in the CellCycleScoring function, which stores in object metadata:
#     * S and G2/M scores
#     * predicted classification of each cell in either G2M, S or G1 phase.
# - subtract ('regress out') this source of heterogeneity from the data.
#    * scaleData(): For each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores.
#    * ScaleData(x, vars.to.regress = c("S.Score", "G2M.Score") ...)

# Solution 2: 
# When analyzing differentiation processes, an alternative workflow should be used
# - Here, regressing out all cell cycle effects can blur the distinction between stem and progenitor cells as well.
# - As an alternative, we suggest regressing out the difference between the G2M and S phase scores.
# - This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle
#   phase amongst proliferating cells (which are often uninteresting), will be regressed out of the data

# In practice:
# - ScaleData(x, vars.to.regress = "CC.Difference" ...)
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #
# Get a list of cell cycle markers (G2/M phase and S phase markers)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Change gene names to uppercase (human genes)
rownames(sce.n) <- toupper(rownames(sce.n))

# Convert sce into Seurat object
data <- as.Seurat(sce.n, counts = "counts", data = "logcounts")

# Perform Seurat normalization
data <- NormalizeData(data)

# Add normalized values from sce object to the Seurat object
data@assays$RNA@data <- logcounts(sce.n)

# Assign cell cycle scores
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Calculate the difference between the G2M and S phase scores
data$CC.Difference <- data$S.Score - data$G2M.Score

# Regress out CC difference
data <- ScaleData(data, vars.to.regress = "CC.Difference", features = rownames(data))

# Save scaled object
saveRDS(data, file = "data_cc")

# Find variable features
data <- FindVariableFeatures(data, selection.method = "vst")

# cell cycle effects strongly mitigated in PCA
data <- RunPCA(data, features = VariableFeatures(data), nfeatures.print = 10)
DimPlot(data)

# PCA on cell cycle genes
data <- RunPCA(data, features = c(s.genes, g2m.genes))
DimPlot(data)


###############################################
### Feature selection & dimension reduction ###
###############################################
# --------------------------------------------------------------- #
# Primary data
# --------------------------------------------------------------- #



# --------------------------------------------------------------- #
# Dell'Orso et al. dataset
# --------------------------------------------------------------- #
# Model per-gene variance (technical & biological variation)
dec.dellOrso <- modelGeneVar(sce.dellOrso.n)

#Visualize the fit
fit1 <- metadata(dec.dellOrso)
plot(fit1$mean, fit1$var, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
curve(fit1$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# Define the highly variable genes (HVGs) and perform dimension reduction
hvgs.dellOrso <- getTopHVGs(dec.dellOrso, prop = 0.1)
sce.dellOrso.n <- runPCA(sce.dellOrso.n, subset_row = hvgs.dellOrso)
plotReducedDim(sce.dellOrso.n, dimred = "PCA", colour_by = "sample")

#TODO: "percentVar" attribute not found.
#percent.var.dellOrso <- attr(reducedDim(sce.dellOrso.n, "percentVar"))
#chosen.elbow.dellOrso <- PCAtools::findElbowPoint(percent.var)
#chosen.elbow.dellOrso
#par(mfrow=c(1,1))
#plot(percent.var.dellOrso, xlab = "PC", ylab = "Variance explained (%)")
#abline(v=chosen.elbow.dellOrso, col = "red")

# Remove PCs corresponding to technical noise
set.seed(123)
denoised.sce.dellOrso <- denoisePCA(sce.dellOrso.n,
                                    technical = dec.dellOrso,
                                    subset.row = hvgs.dellOrso)

# Dimensions of denoised PCA
ncol(reducedDim(denoised.sce.dellOrso))

# Save sce object
saveRDS(denoised.sce.dellOrso, file = "denoised_sce_dellOrso")

# --------------------------------------------------------------- #
# De Micheli et al. dataset
# --------------------------------------------------------------- #
#sce.deMicheli.n <- readRDS("sce_deMicheli_n")

# Model per-gene variance (technical & biological variation)
dec.deMicheli <- modelGeneVar(sce.deMicheli.n)

#Visualize the fit
fit1 <- metadata(dec.deMicheli)
plot(fit1$mean, fit1$var, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
curve(fit1$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# Define the highly variable genes (HVGs) and perform dimension reduction
hvgs.deMicheli <- getTopHVGs(dec.deMicheli, prop = 0.1)
sce.deMicheli.n <- runPCA(sce.deMicheli.n, subset_row = hvgs.deMicheli)
plotReducedDim(sce.deMicheli.n, dimred = "PCA", colour_by = "sample")

#TODO: "percentVar" attribute not found.
#percent.var.deMicheli <- attr(reducedDim(sce.deMicheli.n, "percentVar"))
#chosen.elbow.deMicheli <- PCAtools::findElbowPoint(percent.var)
#chosen.elbow.deMicheli
#par(mfrow=c(1,1))
#plot(percent.var.deMicheli, xlab = "PC", ylab = "Variance explained (%)")
#abline(v=chosen.elbow.deMicheli, col = "red")

# Remove PCs corresponding to technical noise
set.seed(123)
denoised.sce.deMicheli <- denoisePCA(sce.deMicheli.n,
                                    technical = dec.deMicheli,
                                    subset.row = hvgs.deMicheli)

# Dimensions of denoised PCA
ncol(reducedDim(denoised.sce.deMicheli))

# Save sce object
saveRDS(denoised.sce.deMicheli, file = "denoised_sce_deMicheli")


##################
### Clustering ###
##################
# --------------------------------------------------------------- #
# Dell'Orso et al.
# --------------------------------------------------------------- #
# Perform clustering
snng.dellOrso <- buildSNNGraph(sce.dellOrso.n, d = 5)
clust.dellOrso <- igraph::cluster_louvain(snng.dellOrso)$membership
sce.dellOrso.n$cluster <- factor(clust.dellOrso)

set.seed(123)
reducedDim(sce.dellOrso.n, "force") <- igraph::layout_with_fr(snng.dellOrso)
table(clust.dellOrso)

# Plot results
plotReducedDim(sce.dellOrso.n, colour_by = "sample", dimred = "force")
plotReducedDim(sce.dellOrso.n, colour_by = "sample", dimred = "PCA")

sce.dellOrso.n <- runUMAP(sce.dellOrso.n, dimred = "PCA", n_dimred = 5)
plotReducedDim(sce.dellOrso.n, colour_by = "sample", dimred = "UMAP")

# --------------------------------------------------------------- #
# De Micheli et al.
# --------------------------------------------------------------- #
# Perform clustering
snng.deMicheli <- buildSNNGraph(sce.deMicheli.n, d = 5)
clust.deMicheli <- igraph::cluster_louvain(snng.deMicheli)$membership
sce.deMicheli.n$cluster <- factor(clust.deMicheli)

set.seed(123)
reducedDim(sce.deMicheli.n, "force") <- igraph::layout_with_fr(snng.deMicheli)
table(clust.dellOrso)

# Plot results
plotReducedDim(sce.deMicheli.n, colour_by = "sampleID", dimred = "force")
plotReducedDim(sce.deMicheli.n, colour_by = "sampleID", dimred = "PCA")

sce.deMicheli.n <- runUMAP(sce.deMicheli.n, dimred = "PCA", n_dimred = 5)
plotReducedDim(sce.deMicheli.n, colour_by = "sampleID", dimred = "UMAP")




