
# ------------------------------------------------------------------------------
# Aligns gene scores with the CMap rank matrix (keeping only shared genes),
# then ranks gene scores and re-ranks the CMap matrix.
#
# concordance = "negative": rank query ascending (low score = high rank).
# concordance = "positive": rank query descending (high score = high rank).
# ------------------------------------------------------------------------------

align_and_rank <- function(cluster_deg, drug_expr_matrix, concordance) {

  shared_genes           <- intersect(cluster_deg$geneSymbol, rownames(drug_expr_matrix))
  rownames(cluster_deg)  <- cluster_deg$geneSymbol
  cluster_deg            <- cluster_deg[shared_genes, , drop = FALSE]

  cluster_deg$geneRank <- switch(
    concordance,
    "negative" = rank( cluster_deg$score, ties.method = "first"),
    "positive" = rank(-cluster_deg$score, ties.method = "first"),
    stop("concordance must be 'negative' or 'positive'.")
  )

  # Subset to shared genes and rank each instance column
  drug_rank_matrix <- drug_expr_matrix[shared_genes, , drop = FALSE]
  drug_rank_matrix <- as.data.frame(apply(drug_rank_matrix, 2, function(col) rank(-col)))

  return(list(cluster_deg, drug_rank_matrix))
}
