
# ------------------------------------------------------------------------------
# Aggregates per instance p-values into per drug p-values for each cluster
# using 'compute_chisq_pvalues'.
# The computation is parallelized across clusters using 'mclapply()'.
# ------------------------------------------------------------------------------

aggregate_instance_pvalues <- function(drug_info, pval_list, cores) {

  unique_drugs      <- unique(drug_info$cmap_name)
  unique_drug_names <- gsub(".*_", "", unique_drugs)

  cluster_names <- names(pval_list[["pmin"]])

  # Aggregate by cluster
  aggregate_by_cluster <- function(cluster_name) {
    pmin_mat <- as.data.frame(pval_list[["pmin"]][[cluster_name]])
    pmax_mat <- as.data.frame(pval_list[["pmax"]][[cluster_name]])

    n_genes        <- nrow(pmin_mat)
    n_unique_drugs <- length(unique_drugs)

    aggregated_pmin <- as.data.frame(matrix(
      0, nrow = n_genes, ncol = n_unique_drugs,
      dimnames = list(rownames(pmin_mat), unique_drug_names)
    ))
    aggregated_pmax <- aggregated_pmin

    for (drug_idx in seq_along(unique_drugs)) {
      instance_ids <- as.character(
        drug_info[drug_info$cmap_name == unique_drugs[drug_idx], "instance_id"]
      )
      drug_col <- unique_drug_names[drug_idx]

      # pmin
      drug_pmin <- pmin_mat[, instance_ids, drop = FALSE]
      aggregated_pmin[[drug_col]] <- if (length(instance_ids) == 1) {
        drug_pmin[[1]]
      } else {
        apply(drug_pmin, 1, compute_chisq_pvalues)
      }

      # pmax
      drug_pmax <- pmax_mat[, instance_ids, drop = FALSE]
      aggregated_pmax[[drug_col]] <- if (length(instance_ids) == 1) {
        drug_pmax[[1]]
      } else {
        apply(drug_pmax, 1, compute_chisq_pvalues)
      }
    }

    list(pmin = aggregated_pmin, pmax = aggregated_pmax)
  }

  effective_cores <- min(cores, parallel::detectCores() - 1L, length(cluster_names))
  cluster_results <- parallel::mclapply(
    cluster_names, aggregate_by_cluster,
    mc.cores = max(1L, effective_cores)
  )
  names(cluster_results) <- cluster_names

  list(
    aggregated_pmin = setNames(lapply(cluster_results, `[[`, "pmin"), cluster_names),
    aggregated_pmax = setNames(lapply(cluster_results, `[[`, "pmax"), cluster_names)
  )
}
