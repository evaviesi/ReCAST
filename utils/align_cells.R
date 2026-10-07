library(Seurat)
library(SingleR)
library(celldex)

# ------------------------------------------------------------------------------
# This function performs normalization, variable feature selection and integrates
# the input datasets. It optionally regresses out cell cycle effects and performs
# PCA, UMAP and clustering. Cell types are then assigned using SingleR with the
# Human Primary Cell Atlas reference.
# ------------------------------------------------------------------------------

align_cells <- function(SC_list,
                        CellCycle = TRUE,
                        anchor_features = 2000,
                        ndim = 15,
                        npcs = 15,
                        resolution = 0.4,
                        by_CellType = TRUE,
                        seed = 42) {

    set.seed(seed)
    message("Normalizing Data and Finding Variable Features")
    for (i in 1:length(SC_list)) {
    SC_list[[i]] <- NormalizeData(SC_list[[i]], verbose = FALSE)
    SC_list[[i]] <- FindVariableFeatures(SC_list[[i]], selection.method = "vst",
                                         nfeatures = anchor_features, verbose = FALSE)
  }

  message("Integrating Data, these may take some time...")
  SC_anchors <- FindIntegrationAnchors(object.list = SC_list, anchor.features = anchor_features, dims = 1:ndim)
  SC_integrated <- IntegrateData(anchorset = SC_anchors, dims = 1:ndim)
  DefaultAssay(SC_integrated) <- "integrated"

  if (CellCycle) {
    message("Cell Cycle Regression")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    SC_integrated <- CellCycleScoring(SC_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    SC_integrated <- ScaleData(SC_integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC_integrated))
    SC_integrated <- RunPCA(SC_integrated, npcs = npcs, verbose = FALSE)
  }	else {

    # Run the workflow for visualization and clustering
    SC_integrated <- ScaleData(SC_integrated, verbose = FALSE)
    SC_integrated <- RunPCA(SC_integrated, npcs = npcs, verbose = FALSE)
  }

  # Run UMAP and clustering
  SC_integrated <- RunUMAP(SC_integrated, reduction = "pca", dims = 1:ndim)
  SC_integrated <- FindNeighbors(SC_integrated, reduction = "pca", dims = 1:ndim)
  SC_integrated <- FindClusters(SC_integrated, algorithm = 1, resolution = resolution)

  message("Cell Type Annotation")
  if (by_CellType == TRUE) {
    if (class(SC_integrated@assays$RNA) == "Assay") {
      data <- as.matrix(SC_integrated@assays$RNA@data)
      hpca_se <- HumanPrimaryCellAtlasData()
      pred_hpca <- SingleR(test = data, ref = hpca_se, assay.type.test=1, labels = hpca_se$label.main)
      cell_label <- data.frame(row.names = row.names(pred_hpca), celltype=pred_hpca$labels)
    } else {
      datalayers <- Layers(SC_integrated@assays$RNA, search = "data")
      hpca_se <- HumanPrimaryCellAtlasData()
      cell_label <- data.frame(celltype = character(0))
      for (i in 1:length(datalayers)) {
        cat(paste0("Annotating layer: ",  datalayers[i], "\n"))
        data <- as.matrix(SC_integrated@assays$RNA[datalayers[i]])
        pred_hpca <- SingleR(test = data, ref = hpca_se, assay.type.test = 1, labels = hpca_se$label.main)
        cell_label_tmp <- data.frame(row.names = row.names(pred_hpca), celltype = pred_hpca$labels)
        cell_label <- rbind(cell_label, cell_label_tmp)
      }
    }

    if(length(SC_integrated@meta.data$celltype) > 0) {
      SC_integrated@meta.data$celltype <- cell_label$celltype
    } else {
      SC_integrated@meta.data <- cbind(SC_integrated@meta.data, cell_label)
    }

    new_cells <- data.frame()
    for (i in unique(SC_integrated$seurat_clusters)) {
      sub.data <- subset(SC_integrated, seurat_clusters == i)
      temp <- table(sub.data@meta.data$celltype)
      best_cell <- names(which(temp == temp[which.max(temp)]))
      cells_temp <- data.frame(cell_id = row.names(sub.data@meta.data), celltype = best_cell)
      new_cells <- rbind(new_cells, cells_temp)
    }
    cell_meta <- SC_integrated@meta.data
    cell_id <- rownames(cell_meta)
    row.names(new_cells) <- new_cells[, 1]
    new_cells <- new_cells[cell_id, ]
    SC_integrated@meta.data$celltype <- new_cells$celltype
  } else {
    SC_integrated@meta.data$celltype <- paste0("C", as.numeric(SC_integrated@meta.data$seurat_clusters))
  }
  return(SC_integrated)
}
