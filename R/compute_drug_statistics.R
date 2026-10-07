#' @title Compute drug statistics for each cluster
#'
#' @description
#' For each cluster or cell type, computes concordance p-values between
#' differentially expressed genes and drug instance gene profiles in the
#' CMap reference, then aggregates to the drug level via a one-sided KS test.
#'
#' The computation is parallelized across clusters using 'mclapply()'.
#'
#' @param cluster_degs Named list with a data.frame of differentially expressed
#' genes per cluster.
#' @param drug_ref_profiles Output of 'generate_reference' function.
#' @param repurposing_unit  Either 'drug' (default) or 'treatment'.
#' @param CEG_threshold Significance threshold for concordantly expressed genes
#' (default 0.05).
#' @param concordance Direction of concordance: either 'negative' or 'positive'.
#' @param drug_type One of 'FDA', 'compounds', or 'ALL'.
#' @param imp_genes Logical; whether to return per drug gene-level p-values.
#' @param cores Number of cores used to run the process in parallel.
#'
#' @return When imp_genes = FALSE: a named list of drug p-value data.frames,
#' one per cluster.  When imp_genes = TRUE: a list with 'drug' and 'pvalue'.
#' @export

compute_drug_statistics <- function(cluster_degs       = NULL,
                                    drug_ref_profiles  = NULL,
                                    repurposing_unit   = "drug",
                                    CEG_threshold      = 0.05,
                                    concordance        = "negative",
                                    drug_type          = "FDA",
                                    imp_genes          = FALSE,
                                    cores              = 1) {

  # Input validation
  if (!(drug_type %in% c("FDA", "compounds", "ALL"))) {
    stop("drug_type must be 'FDA', 'compounds', or 'ALL'.")
  }
  if (!(repurposing_unit %in% c("drug", "treatment"))) {
    stop("repurposing_unit must be 'drug' or 'treatment'.")
  }

  # Filter reference profiles
  if (drug_type != "ALL") {
    drug_names <- gsub("_.*", "", drug_ref_profiles$drug_info$cmap_name)
    if (drug_type == "FDA") {
      keep <- drug_names %in% FDA_drugs
    } else {
      keep <- !(drug_names %in% FDA_drugs)
    }
    drug_ref_profiles$drug_info        <- drug_ref_profiles$drug_info[keep, , drop = FALSE]
    drug_ref_profiles$drug_expr_matrix <- drug_ref_profiles$drug_expr_matrix[
      , as.character(drug_ref_profiles$drug_info$instance_id), drop = FALSE
    ]
  }

  n_clusters      <- length(cluster_degs)
  effective_cores <- max(1, min(cores, parallel::detectCores() - 1, n_clusters))
  message("Running on ", effective_cores, " core(s) (max available cores: ",
          parallel::detectCores(), ", max parallelizable clusters: ", n_clusters, ")")

  drug_info <- drug_ref_profiles$drug_info
  drug_info$drug <- if (repurposing_unit == "drug") drug_info$cmap_name else drug_info$treatment

  # ============================================================
  # Step 1 — Per gene beta p-values for each cluster
  # ============================================================
  compute_gene_pvals_for_cluster <- function(cluster_name) {
    cluster_deg <- data.frame(
      geneSymbol = rownames(cluster_degs[[cluster_name]]),
      score      = cluster_degs[[cluster_name]]$score,
      stringsAsFactors = FALSE
    )

    aligned_data     <- align_and_rank(cluster_deg, drug_ref_profiles$drug_expr_matrix, concordance)
    aligned_deg      <- aligned_data[[1]]
    aligned_rank_mat <- aligned_data[[2]]

    return(list(
      pmin = compute_beta_pvalues("min", aligned_rank_mat, aligned_deg),
      pmax = compute_beta_pvalues("max", aligned_rank_mat, aligned_deg)
    ))
  }

  gene_pval_results <- parallel::mclapply(
    names(cluster_degs), compute_gene_pvals_for_cluster,
    mc.cores = effective_cores
  )
  names(gene_pval_results) <- names(cluster_degs)

  pmin_per_cluster <- setNames(lapply(gene_pval_results, `[[`, "pmin"), names(cluster_degs))
  pmax_per_cluster <- setNames(lapply(gene_pval_results, `[[`, "pmax"), names(cluster_degs))

  # ============================================================
  # Step 2 — Per drug concordance p-values for each cluster
  # ============================================================
  compute_drug_pvals_for_cluster <- function(cluster_name) {

    pmin_mat <- pmin_per_cluster[[cluster_name]]
    pmax_mat <- pmax_per_cluster[[cluster_name]]

    # Take the most significant orientation at each gene
    min_gene_pval  <- pmin(pmin_mat, pmax_mat)
    z_scores       <- qnorm(min_gene_pval, lower.tail = FALSE)

    # CEG score: sum of Z-scores above the significance threshold
    z_threshold  <- qnorm(CEG_threshold, lower.tail = FALSE)
    ceg_scores   <- vapply(seq_len(ncol(z_scores)), function(col_idx) {
      col_z <- z_scores[, col_idx]
      sum(col_z[col_z >= z_threshold])
    }, numeric(1))
    names(ceg_scores) <- colnames(z_scores)

    drug_pvals <- compute_ks_pvalues(ceg_scores, drug_info)
    drug_pvals$p_val <- pmin(1, drug_pvals$p_val + 1e-6)
    drug_pvals <- drug_pvals[order(drug_pvals$p_val), ]
    drug_names <- rownames(drug_pvals)

    drug_pvals$drug_name <- gsub("_BRD-.*", "", drug_names)
    drug_pvals$drug_id   <- gsub(".*_",     "", drug_names)
    rownames(drug_pvals) <- NULL

    drug_pvals     <- drug_pvals[, c("drug_name", "drug_id", "p_val")]
    drug_pvals$fdr <- p.adjust(drug_pvals$p_val, method = "fdr")

    return(drug_pvals)
  }


  drug_pval_results <- parallel::mclapply(
    names(cluster_degs), compute_drug_pvals_for_cluster,
    mc.cores = effective_cores
  )

  names(drug_pval_results) <- names(cluster_degs)
  result <- list(drug = drug_pval_results)

  # ============================================================
  # Step 3 (optional) — Aggregate per instance p-values for
  # important gene identification
  # ============================================================
  if (imp_genes) {
    pval_list <- list(pmin = pmin_per_cluster, pmax = pmax_per_cluster)

    if (repurposing_unit == "drug") {
      result$pvalue <- aggregate_instance_pvalues(drug_info, pval_list, cores)
    } else {
      result$pvalue <- pval_list
    }
  }

  if (imp_genes) return(result) else return(result$drug)
}
