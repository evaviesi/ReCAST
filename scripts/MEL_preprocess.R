# Load packages and utilities --------------------------------------------------

library(Seurat)
library(qs2)
library(rlang)
library(rhdf5)
source("utils/align_cells.R")


# Set parameters ---------------------------------------------------------------

set.seed(42)
data_path <- "raw_data"
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
out_path <- "out_data"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


# Load data --------------------------------------------------------------------

## Download and load patient scRNA-seq data (GEO: GSE189889)
h5file <- file.path(data_path, "GSE189889_acral_2101.debatched.iscva.h5")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE189889&format=file&file=GSE189889%5Facral%5F2101%2Edebatched%2Eiscva%2Eh5",
              destfile = h5file,
              mode = "wb")
h5ls(h5file)

## Read components
barcodes   <- h5read(h5file, "matrix/barcodes")
gene_names <- h5read(h5file, "matrix/gene_names")
data       <- h5read(h5file, "matrix/data")
indices    <- h5read(h5file, "matrix/indices")
indptr     <- h5read(h5file, "matrix/indptr")
shape      <- h5read(h5file, "matrix/shape")

## Construct sparse matrix
melanoma_counts <- new("dgCMatrix",
                       Dim = as.integer(shape),
                       Dimnames = list(gene_names, barcodes),
                       x = as.numeric(data),
                       i = as.integer(indices),
                       p = as.integer(indptr))

## Subset data retaining only primary tumor cells (remove metastatic samples)
sample_to_remove <- c("Acral-01", "MOF-AM5", "MOF-AM7-RNA","MOF-AM8-Node-RNA")
melanoma_metadata <- h5read(h5file, "artifacts")
prim_tumor_id <- melanoma_metadata$all$covs[!(melanoma_metadata$all$covs$sample %in% sample_to_remove),"id"]
metadata <- melanoma_metadata$all$covs[!(melanoma_metadata$all$covs$sample %in% sample_to_remove), ]
rownames(metadata) <- metadata$id
melanoma_counts_prim <- as.matrix(melanoma_counts)[,prim_tumor_id]
dim(melanoma_counts_prim) # [1] 33538 23469

## Create a Seurat object for each sample
samples <- unique(metadata$sample)
SC_melanoma_list <- list()
for (sample in samples) {
  cat(paste0("Load sample: ", sample, "\n"))
  metadata_tmp <- metadata[metadata$sample == sample, ]
  counts_tmp <- melanoma_counts_prim[ , rownames(metadata_tmp)]
  SC_melanoma_list[sample] <- CreateSeuratObject(counts       = counts_tmp,
                                                 project      = "MEL",
                                                 min.cells    = 3,
                                                 min.features = 200,
                                                 meta.data    = metadata_tmp)
  SC_melanoma_list[[sample]]$orig.ident <- sample
}


## Download and load control scRNA-seq data (GEO: GSE183047)
file_path <- file.path(data_path, "GSE183047_RAW.tar")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183047&format=file",
              destfile = file_path,
              mode = "wb")
file_names <- untar(file_path, list = TRUE)

ctr_samples <- c("GSM5550119", "GSM5550120", "GSM5550121", "GSM5550122",
                 "GSM5550123", "GSM5550124", "GSM5550125", "GSM5550126",
                 "GSM5550127", "GSM5550128")

for (ctr_sample in ctr_samples) {
  files_to_extract <- file_names[grepl(ctr_sample, file_names)]
  # extract files
  untar(file_path, files = files_to_extract, exdir = paste0(data_path, "/", ctr_sample))
  # rename files
  extracted_files <- list.files(paste0(data_path, "/", ctr_sample), full.names = TRUE)
  file.rename(extracted_files, file.path(dirname(extracted_files), sub("^.*_", "", basename(extracted_files))))
}

## Create a Seurat object for each sample
SC_control_list <- list()
for (file in ctr_samples) {
  cat(paste0("Load sample: ", file, "\n"))
  counts_tmp <- Read10X(data.dir = file.path(data_path, file))
  metadata_tmp <- data.frame(row.names = colnames(counts_tmp),
                             id        = colnames(counts_tmp),
                             sample    = file)
  SC_control_list[file] <- CreateSeuratObject(counts = counts_tmp,
                                              project = "MEL",
                                              min.cells = 3,
                                              min.features = 200,
                                              meta.data = metadata_tmp)
  SC_control_list[[file]]$orig.ident <- file
}


# QC and filtering -------------------------------------------------------------

## Merge Seurat objects
SC_list <- c(SC_melanoma_list, SC_control_list)
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
sample[which(sample == "GSM5550119")] <- "Normal1"
sample[which(sample == "GSM5550120")] <- "Normal2"
sample[which(sample == "GSM5550121")] <- "Normal3"
sample[which(sample == "GSM5550122")] <- "Normal4"
sample[which(sample == "GSM5550123")] <- "Normal5"
sample[which(sample == "GSM5550124")] <- "Normal6"
sample[which(sample == "GSM5550125")] <- "Normal7"
sample[which(sample == "GSM5550126")] <- "Normal8"
sample[which(sample == "GSM5550127")] <- "Normal9"
sample[which(sample == "GSM5550128")] <- "Normal10"
SC_merge@meta.data$sample <- sample

## Filtering by 95% quantile mitochondrial percentage
mito_thresh <- tapply(SC_merge$percent.mito, SC_merge$sample, quantile, 0.95)
SC_merge <- subset(SC_merge, subset =
                     (sample == "Normal1" & percent.mito < 15) |
                     (sample == "Normal2" & percent.mito < 15) |
                     (sample == "Normal3" & percent.mito < 15) |
                     (sample == "Normal4" & percent.mito < 15) |
                     (sample == "Normal5" & percent.mito < 15) |
                     (sample == "Normal6" & percent.mito < 15) |
                     (sample == "Normal7" & percent.mito < 15) |
                     (sample == "Normal8" & percent.mito < 15) |
                     (sample == "Normal9" & percent.mito < 15) |
                     (sample == "Normal10" & percent.mito < 15)|
                     (sample == "MOF-AM2" & percent.mito < mito_thresh["MOF-AM2"])         |
                     (sample == "MOF-AM3" & percent.mito < mito_thresh["MOF-AM3"])         |
                     (sample == "MOF-AM4" & percent.mito < mito_thresh["MOF-AM4"])         |
                     (sample == "MOF-AM6-RNA" & percent.mito < mito_thresh["MOF-AM6-RNA"]) |
                     (sample == "MOF-AM8-Toe-RNA" & percent.mito < mito_thresh["MOF-AM8-Toe-RNA"]))

## Filtering by 99% quantile ribosomal percentage
ribo_thresh <- tapply(SC_merge$percent.ribo, SC_merge$sample, quantile, 0.99)
SC_merge <- subset(SC_merge, subset =
                     (sample == "Normal1" & percent.ribo < ribo_thresh["Normal1"]) |
                     (sample == "Normal2" & percent.ribo < ribo_thresh["Normal2"]) |
                     (sample == "Normal3" & percent.ribo < ribo_thresh["Normal3"]) |
                     (sample == "Normal4" & percent.ribo < ribo_thresh["Normal4"]) |
                     (sample == "Normal5" & percent.ribo < ribo_thresh["Normal5"]) |
                     (sample == "Normal6" & percent.ribo < ribo_thresh["Normal6"]) |
                     (sample == "Normal7" & percent.ribo < ribo_thresh["Normal7"]) |
                     (sample == "Normal8" & percent.ribo < ribo_thresh["Normal8"]) |
                     (sample == "Normal9" & percent.ribo < ribo_thresh["Normal9"]) |
                     (sample == "Normal10" & percent.ribo < ribo_thresh["Normal10"])       |
                     (sample == "MOF-AM2" & percent.ribo < ribo_thresh["MOF-AM2"])         |
                     (sample == "MOF-AM3" & percent.ribo < ribo_thresh["MOF-AM3"])         |
                     (sample == "MOF-AM4" & percent.ribo < ribo_thresh["MOF-AM4"])         |
                     (sample == "MOF-AM6-RNA" & percent.ribo < ribo_thresh["MOF-AM6-RNA"]) |
                     (sample == "MOF-AM8-Toe-RNA" & percent.ribo < ribo_thresh["MOF-AM8-Toe-RNA"]))

## Filtering by the number of detected features and total counts
SC_merge <- subset(SC_merge, nCount_RNA > 500 & nFeature_RNA > 200 & nCount_RNA < 50000 & nFeature_RNA < 7000)
VlnPlot(SC_merge, features = "percent.mito", group.by = "sample")
VlnPlot(SC_merge, features = "percent.ribo", group.by = "sample")
FeatureScatter(SC_merge, "nCount_RNA", "nFeature_RNA", group.by = "sample")

## Remove low-quality control samples
table(SC_merge$sample)
SC_merge <- subset(SC_merge, !(sample %in%  c("Normal2", "Normal3", "Normal4", "Normal5", "Normal6")))


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
