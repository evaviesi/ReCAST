
# ------------------------------------------------------------------------------
# Computes a p-value for each drug via a one-sided KS test:
# the drug's own instances vs. all remaining instances.
# ------------------------------------------------------------------------------

compute_ks_pvalues <- function(ceg_scores, drug_info) {

  unique_drug_names <- unique(drug_info$drug)

  suppressWarnings(
    pval_per_drug <- vapply(unique_drug_names, function(drug_id) {
      instance_ids      <- as.character(drug_info[drug_info$drug == drug_id, "instance_id"])
      instance_scores   <- ceg_scores[instance_ids]
      background_scores <- ceg_scores[setdiff(names(ceg_scores), instance_ids)]
      ks.test(instance_scores, background_scores, alternative = "less")$p.value
    }, numeric(1))
  )

  ks_df <- data.frame(
    p_val     = pval_per_drug,
    drug      = sapply(unique_drug_names, function(x) strsplit(x, "_")[[1]][1]),
    row.names = unique_drug_names,
    stringsAsFactors = FALSE
  )
  return(ks_df)
}

