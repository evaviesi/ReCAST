# Load packages and utilities --------------------------------------------------

library(Seurat)
library(qs2)
library(rlang)
source("utils/align_cells.R")


# Set parameters ---------------------------------------------------------------

set.seed(42)
data_path <- "raw_data"
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
out_path <- "out_data"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


# Load data --------------------------------------------------------------------

## Download and load patient and control scRNA-seq data (GEO: GSE161529)
file_path <- file.path(data_path, "GSE161529_RAW.tar")
options(timeout = 3600)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161529&format=file",
              destfile = file_path,
              mode = "wb")
file_names <- untar(file_path, list = TRUE)

## Download features
feature_path <- file.path(data_path, "features.tsv.gz")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161529&format=file&file=GSE161529%5Ffeatures%2Etsv%2Egz",
              destfile = feature_path,
              mode = "wb")

samples <- c("GSM4909253", "GSM4909254", "GSM4909261", "GSM4909263",
             "GSM4909281", "GSM4909282", "GSM4909283", "GSM4909284")

for (sample in samples) {
  files_to_extract <- file_names[grepl(sample, file_names)]
  # extract files
  untar(file_path, files = files_to_extract, exdir = paste0(data_path, "/", sample))
  # rename files
  extracted_files <- list.files(paste0(data_path, "/", sample), full.names = TRUE)
  file.rename(extracted_files, file.path(dirname(extracted_files), sub("^.*-", "", basename(extracted_files))))
  # add features
  file.copy(from = feature_path, to = file.path(data_path, sample, basename(feature_path)), overwrite = FALSE)
}

## Create a Seurat object for each sample
SC_list <- list()
for (file in samples){
  message(paste0("Loading sample: ", file, "\n"))
  counts_tmp <- Read10X(data.dir = file.path(data_path, file))
  metadata_tmp <- data.frame(row.names = colnames(counts_tmp),
                             cell = colnames(counts_tmp),
                             sample = file)
  SC_list[file] <- CreateSeuratObject(counts       = counts_tmp,
                                      project      = "TNBC",
                                      min.cells    = 3,
                                      min.features = 200,
                                      meta.data    = metadata_tmp)
  SC_list[[file]]$orig.ident <- file
}


# QC and filtering -------------------------------------------------------------

## Merge Seurat objects
SC_merge <- merge(SC_list[[1]], SC_list[-1])

## Compute and visualize QC metrics
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^MT-", col.name = "percent.mito")
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^RP[SL]", col.name = "percent.ribo")
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^HB[^(P)]", col.name = "percent.hemo")

VlnPlot(SC_merge, features = "percent.mito", group.by = "sample")
VlnPlot(SC_merge, features = "percent.ribo", group.by = "sample")
VlnPlot(SC_merge, features = "percent.hemo", group.by = "sample")
FeatureScatter(SC_merge, "nCount_RNA", "nFeature_RNA", group.by = "sample")

## Change sample names
sample <- SC_merge@meta.data$sample
sample[which(sample == "GSM4909253")] <- "Normal1"
sample[which(sample == "GSM4909254")] <- "Normal2"
sample[which(sample == "GSM4909261")] <- "Normal3"
sample[which(sample == "GSM4909263")] <- "Normal4"
sample[which(sample == "GSM4909281")] <- "TNBC1"
sample[which(sample == "GSM4909282")] <- "TNBC2"
sample[which(sample == "GSM4909283")] <- "TNBC3"
sample[which(sample == "GSM4909284")] <- "TNBC4"
SC_merge@meta.data$sample <- sample

## Filtering by 95% quantile mitochondrial percentage
mito_thresh <- tapply(SC_merge$percent.mito, SC_merge$sample, quantile, 0.95)
SC_merge <- subset(SC_merge, subset =
                     (sample == "Normal1" & percent.mito < mito_thresh["Normal1"]) |
                     (sample == "Normal2" & percent.mito < mito_thresh["Normal2"]) |
                     (sample == "Normal3" & percent.mito < mito_thresh["Normal3"]) |
                     (sample == "Normal4" & percent.mito < mito_thresh["Normal4"]) |
                     (sample == "TNBC1" & percent.mito < mito_thresh["TNBC1"])     |
                     (sample == "TNBC2" & percent.mito < mito_thresh["TNBC2"])     |
                     (sample == "TNBC3" & percent.mito < mito_thresh["TNBC3"])     |
                     (sample == "TNBC4" & percent.mito < mito_thresh["TNBC4"]))

## Filtering by the number of detected features and total counts
SC_merge <- subset(SC_merge, nCount_RNA > 500 & nFeature_RNA > 200 & nCount_RNA < 100000 & nFeature_RNA < 7500)
VlnPlot(SC_merge, features = "percent.mito", group.by = "sample")
FeatureScatter(SC_merge, "nCount_RNA", "nFeature_RNA", group.by = "sample")


# Align and annotate cells -----------------------------------------------------

SC_list <- SplitObject(SC_merge, split.by = "sample")
SC_integrated <- align_cells(SC_list = SC_list,
                             CellCycle = TRUE,
                             anchor_features = 2000,
                             ndim = 15,
                             npcs = 15,
                             resolution = 0.4,
                             by_CellType = TRUE)


# Plot and save results --------------------------------------------------------

DimPlot(SC_integrated,
        reduction = "umap",
        label = FALSE,
        group.by = "celltype",
        split.by = "sample")

qs_save(SC_integrated, file = file.path(out_path, "SC_obj.qs2"))
