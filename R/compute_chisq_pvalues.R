
# ------------------------------------------------------------------------------
# Fisher's combined probability test: returns a single combined p-value
# from a vector of independent p-values.
# ------------------------------------------------------------------------------

compute_chisq_pvalues <- function(p_values) {
  valid_pvals <- (p_values > 0) & (p_values <= 1)
  n_valid     <- sum(valid_pvals)

  if (n_valid < 2) {
    warning("compute_chisq_pvalues: need at least two valid p-values; returning NA.")
    return(NA_real_)
  }

  if (n_valid < length(p_values)) {
    warning("compute_chisq_pvalues: some p-values were out of (0, 1] and were omitted.")
  }

  chisq_stat  <- -2 * sum(log(p_values[valid_pvals]))
  df          <- 2 * n_valid
  res_stat    <- pchisq(chisq_stat, df, lower.tail = FALSE)
  return(res_stat)
}
