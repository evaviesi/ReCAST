
# ------------------------------------------------------------------------------
# Computes p-values per gene and per instance using the order-statistic
# beta distribution.
#
# rank_order = "min": tests whether the minimum of (query rank, CMap rank)
# is smaller than expected (detects genes ranked low in both).
# rank_order = "max": tests whether the maximum of (query rank, CMap rank)
# is smaller than expected (detects genes ranked high in both).
# ------------------------------------------------------------------------------

compute_beta_pvalues <- function(rank_order, aligned_rank_mat, aligned_deg) {

  n_genes <- nrow(aligned_rank_mat)
  n_sig   <- ncol(aligned_rank_mat)

  deg_rank_matrix <- matrix(
    aligned_deg$geneRank,
    nrow = n_genes,
    ncol = n_sig
  )

  if (rank_order == "min") {
    order_statistics <- pmin(as.matrix(aligned_rank_mat), deg_rank_matrix)
    p_values <- 1 - pbeta((order_statistics - 1) / n_genes, 1, 2, lower.tail = TRUE)
  } else if (rank_order == "max") {
    order_statistics <- pmax(as.matrix(aligned_rank_mat), deg_rank_matrix)
    p_values <- pbeta(order_statistics / n_genes, 2, 1, lower.tail = TRUE)
  } else {
    stop("rank_order must be 'min' or 'max'.")
  }
  p_values <- p_values + 1e-6
  p_values[p_values > 1] <- 1

  return(p_values)
}
