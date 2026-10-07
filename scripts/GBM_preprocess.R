# Load packages ----------------------------------------------------------------

library(Seurat)
library(qs2)
library(rlang)


# Set parameters ---------------------------------------------------------------

set.seed(42)
data_path <- "raw_data"
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
out_path <- "out_data"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


# Load data --------------------------------------------------------------------

# Load Seurat objects SNU24.rds, SNU25.rds and SNU57normal.rds downloaded
# from the SNUH-GliomaAtlas dataset: https://gbmvisium.snu.ac.kr/ --------------

pat1 <- readRDS(file.path(data_path, "SNU24.rds"))
pat2 <- readRDS(file.path(data_path, "SNU25.rds"))
ctr  <- readRDS(file.path(data_path, "SNU57normal.rds"))


# QC and filtering -------------------------------------------------------------

## Merge Seurat objects
SC_merge <- merge(pat1, c(pat2, ctr))

## Compute and visualize QC metrics
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^MT-", col.name = "percent.mito")
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^RP[SL]", col.name = "percent.ribo")
SC_merge <- PercentageFeatureSet(SC_merge, pattern = "^HB[^(P)]", col.name = "percent.hemo")

VlnPlot(SC_merge, features = "percent.mito", group.by = "sample")
VlnPlot(SC_merge, features = "percent.ribo", group.by = "sample")
VlnPlot(SC_merge, features = "percent.hemo", group.by = "sample")
FeatureScatter(SC_merge, "nCount_RNA", "nFeature_RNA", group.by = "sample")

## Filtering by mitochondrial percentage
SC_merge <- subset(SC_merge, subset =
                     (sample == "SNU24" & percent.mito < 10) |
                     (sample == "SNU25" & percent.mito < 15) |
                     (sample == "SNU57normal" & percent.mito < 10))

## Filtering by 99% quantile ribosomal and hemo percentage
ribo_thresh <- tapply(SC_merge$percent.ribo, SC_merge$sample, quantile, 0.99)
hemo_thresh <- tapply(SC_merge$percent.hemo, SC_merge$sample, quantile, 0.99)
SC_merge <- subset(SC_merge, subset =
                     (sample == "SNU24" & percent.ribo < ribo_thresh["SNU24"]) & percent.hemo < hemo_thresh["SNU24"] |
                     (sample == "SNU25" & percent.ribo < ribo_thresh["SNU25"]) & percent.hemo < hemo_thresh["SNU25"] |
                     (sample == "SNU57normal" & percent.ribo < ribo_thresh["SNU57normal"] & percent.hemo < hemo_thresh["SNU57normal"]))

## Filtering by the number of detected features and total counts
SC_merge <- subset(SC_merge, nCount_RNA > 500 & nFeature_RNA > 200 & nCount_RNA < 100000 & nFeature_RNA < 10000)
VlnPlot(SC_merge, features = "percent.mito", group.by = "sample")
VlnPlot(SC_merge, c("percent.ribo"), group.by = "sample")
VlnPlot(SC_merge, c("percent.hemo"), group.by = "sample")
FeatureScatter(SC_merge, "nCount_RNA", "nFeature_RNA", group.by = "sample")


# Preprocessing ----------------------------------------------------------------

SC_merge <- NormalizeData(SC_merge, scale.factor = 1e6)
SC_merge <- FindVariableFeatures(SC_merge, selection.method = "vst", nfeatures = 2000)
SC_merge <- ScaleData(SC_merge, features = rownames(SC_merge))
SC_merge <- RunPCA(SC_merge, features = VariableFeatures(SC_merge), npcs = 12)
ElbowPlot(SC_merge, ndims = 50)
SC_merge <- RunUMAP(SC_merge, reduction = "pca", dims = 1:12)

## Standardize the column names in metadata
annotation <- "celltype"
SC_merge@meta.data[[annotation]] <- SC_merge@meta.data$finalCelltype_15Sept
SC_merge@meta.data$disease <- SC_merge@meta.data$histology

## Remove the "Doublet" cell type and cell types with < 3 diseased or control cells
SC_merge_subset <- subset(x = SC_merge, subset = !!sym(annotation) != "Doublet")
meta_data <- SC_merge_subset@meta.data
counts <- table(meta_data[[annotation]], meta_data$disease)
remain_cell_types <- rownames(counts)[apply(counts, 1, function(x) all(x >= 3))]
SC_merge_subset <- subset(SC_merge_subset, !!sym(annotation) %in% remain_cell_types)


# Plot and save results --------------------------------------------------------

DimPlot(SC_merge_subset,
        reduction = "umap",
        label = FALSE,
        group.by = "celltype",
        split.by = "disease")

qs_save(SC_merge_subset, file = file.path(out_path, "SC_obj.qs2"))
