#' @title Compute drug-2-cell matrix
#'
#' @description
#' Scores each drug for each cell using the Seurat background-corrected scoring
#' approach: raw drug-target expression scores are adjusted by the mean score
#' of a random control gene set drawn from the same expression-level bins.
#'
#' @param SC_obj A Seurat object containing normalized RNA counts.
#' @param common_genes Character vector of candidate drug-target genes.
#' @param target_ppi_df A data.frame with columns 'drug_name' and 'gene_names'.
#' @param n_bins Number of expression bins for background sampling (default 25).
#' @param ctrl_size Maximum number of control genes drawn per bin (default 50).
#' @param seed Random seed for reproducibility (default 42).
#' @param out_path Output directory path.
#' @param out_name Output file name.
#'
#' @return Drug-2-cell matrix (drugs x clusters/cell types).
#' @export
compute_d2c_matrix <- function(SC_obj,
                               common_genes,
                               target_ppi_df,
                               n_bins    = 25,
                               ctrl_size = 50,
                               seed      = 42,
                               out_path  = "./",
                               out_name  = "d2c_matrix") {

  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # Get expression matrix for drug targets
  SC_obj_subset    <- subset(x = SC_obj, features = common_genes)
  cell_genes_mat <- t(as.matrix(SC_obj_subset@assays$RNA$data))

  drug_names <- target_ppi_df$drug_name

  # Build binary weight matrix
   target_mask_mat <- vapply(drug_names, function(drug_name) {
    target_genes <- strsplit(
      target_ppi_df[target_ppi_df$drug_name == drug_name, "gene_names"],
      ";"
    )[[1]]
    common_genes %in% target_genes
  }, logical(length(common_genes)))
  rownames(target_mask_mat) <- common_genes

  # Normalize each drug column by the number of its targets
  drug_weight_mat <- sweep(target_mask_mat, 2, colSums(target_mask_mat) + 1e-6, "/")

  # Get drug signature per cell
  orig_drug_scores <- cell_genes_mat %*% drug_weight_mat

  # Background correction
  gene_mean_expression <- colMeans(cell_genes_mat)
  n_items_per_bin      <- round(length(gene_mean_expression) / (n_bins - 1))
  gene_bin_assignment  <- rank(gene_mean_expression, ties.method = "min") - 1
  gene_bin_assignment  <- gene_bin_assignment %/% n_items_per_bin

  set.seed(seed)

  # For each bin, randomly select up to ctrl_size genes as controls
  sample_control_genes <- function(bin_id) {
    control_mask <- (gene_bin_assignment == bin_id)
    bin_indices <- which(control_mask)
    sampled_indices  <- sample(bin_indices)
    if (length(sampled_indices) > ctrl_size) {
      control_mask[] <- FALSE
      control_mask[sampled_indices[seq_len(ctrl_size)]] <- TRUE
    }
    control_mask
  }

  unique_bins <- unique(gene_bin_assignment)
  control_gene_mask_mat <- vapply(
    unique_bins, sample_control_genes,
    logical(length(common_genes))
  )
  colnames(control_gene_mask_mat) <- paste0("bin", unique_bins)

  # Normalize control columns
  control_weight_mat <- sweep(control_gene_mask_mat, 2, colSums(control_gene_mask_mat) + 1e-6, "/")

  # Background profiles per cell
  control_profiles <- cell_genes_mat %*% control_weight_mat

  # Map each drug to the expression bins of its own targets
  drug_bin_mat <- vapply(drug_names, function(drug_name) {
    target_mask      <- target_mask_mat[, drug_name]
    bin_ids          <- unique(gene_bin_assignment[target_mask])
    colnames(control_weight_mat) %in% paste0("bin", bin_ids)
  }, logical(ncol(control_weight_mat)))
  rownames(drug_bin_mat) <- colnames(control_weight_mat)

  # Normalize drug bin columns
  drug_bin_weight_mat <- sweep(drug_bin_mat, 2, colSums(drug_bin_mat) + 1e-6, "/")

  # Background drug score per cell
  background_drug_scores <- control_profiles %*% drug_bin_weight_mat

  # Background-corrected d2c matrix
  d2c_matrix <- t(orig_drug_scores - background_drug_scores)

  saveRDS(d2c_matrix, file.path(out_path, paste0(out_name, ".rds")))
  return(d2c_matrix)
}
