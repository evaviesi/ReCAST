
# ------------------------------------------------------------------------------
# This function computes differentially expressed genes (DEGs) between case and
# control samples for each cluster/cell type using the Seurat Wilcoxon test.
# ------------------------------------------------------------------------------

get_de_genes <- function(SC_obj,
                         annotation = "celltype",
                         seed = 42) {

  set.seed(seed)
  cluster_degs  <- list()
  cluster_names <- c()

  # Subset for each cluster/cell type
  for (i in unique(SC_obj@meta.data[[annotation]])) {

    Idents(SC_obj) <- annotation

    if (annotation == "celltype") {
      cluster_cells <- subset(SC_obj, celltype == i)
    } else {
      cluster_cells <- subset(SC_obj, !!rlang::sym(annotation) == i)
    }

    Idents(cluster_cells) <- "type"

    # ---- Seurat ----
    degs <- Seurat::FindMarkers(cluster_cells, ident.1 = "Case", ident.2 = "Control")
    degs_for_drug <- data.frame(row.names   = row.names(degs),
                                  score     = degs$avg_log2FC,
                                  adj_p_val = degs$p_val_adj,
                                  p_val     = degs$p_val)
    cluster_degs[[i]] <- degs_for_drug
    cluster_names     <- c(cluster_names, i)

  }

  names(cluster_degs) <- cluster_names
  return(cluster_degs)
}
