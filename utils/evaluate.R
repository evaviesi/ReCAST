library(ROCR)
library(PRROC)

# ------------------------------------------------------------------------------
# This function computes ROC and PRC curves and returns the corresponding
# AUC and AUPRC values.
# ------------------------------------------------------------------------------

evaluate <- function(drug_score_df, positive_drugs) {

  drug_score_df$label <- ifelse(tolower(rownames(drug_score_df)) %in% tolower(positive_drugs), 1, 0)

  pred      <- prediction(predictions = drug_score_df$drug_score, labels = drug_score_df$label)
  auc_perf  <- performance(pred, measure = "auc")
  auc_value <- auc_perf@y.values[[1]]

  pos_scores <- drug_score_df$drug_score[drug_score_df$label == 1]
  neg_scores <- drug_score_df$drug_score[drug_score_df$label == 0]

  pr_result  <- pr.curve(scores.class0 = pos_scores,
                               scores.class1 = neg_scores,
                               curve = FALSE)
  aupr_value <- pr_result$auc.integral

  return(c(auc_value, aupr_value))
}
