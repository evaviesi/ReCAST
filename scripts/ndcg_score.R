# Libraries and Packages --------------------------------------------------
library(ReCAST)
library(ggplot2)
library(patchwork)
# Function ----------------------------------------------------------------

# DCG function
dcg <- function(labels, scores, k = length(labels)) {

  ord <- order(scores, decreasing = TRUE)
  rel <- labels[ord][1:k]
  sum((2^rel - 1) / log2(seq_along(rel) + 1))
}

# NDCG function
ndcg <- function(labels, scores, k = length(labels)) {

  actual <- dcg(labels, scores, k)
  ideal_order <- order(labels, decreasing = TRUE)
  ideal_rel <- labels[ideal_order][1:k]
  ideal <- sum((2^ideal_rel - 1) /
                 log2(seq_along(ideal_rel) + 1))
  actual / ideal
}

# Filter positive drugs for disease of interest
load_approved <- function(
    disease   = NULL
) {
  rownames(ReCAST::positive_drugs)[ReCAST::positive_drugs[[disease]] == 1]
}

# Filter score data by common perturbagens, and add labels
load_scores <- function(
    score_path     = "data",
    disease        = NULL,
    approved_drugs = NULL,
    type           = "overall"
) {
  
  # Define score path
  data_path   <- file.path(score_path, disease)
  if (type == "overall") {
    score_files <- list(
      ASGARD     = file.path("ASGARD", paste0(disease, "_drug_scores.csv")),
      scDrugLink = file.path("scDrugLink", paste0(disease, "_drug_scores.csv")),
      ReCAST     = file.path("output", "neg_drug_scores.csv")
    )
  } else {
    score_files <- list(
      ASGARD     = file.path("ASGARD", paste0(disease, "_individual_cell_type_drug_scores.csv")),
      scDrugLink = file.path("scDrugLink", paste0(disease, "_individual_cell_type_drug_scores.csv")),
      ReCAST     = file.path("output", "neg_cluster_scores.csv")
    )
  }
  
  # Load and normalize overall results
  score_list <- lapply(score_files, function(rel_path) {
    score_df <- read.csv(file.path(data_path, rel_path), row.names = 1)
    if (type == "overall") {
      colnames(score_df)[1:3] <- c("drug_score", "p_val", "fdr")
    }
    rownames(score_df) <- tolower(gsub("\\.", "-", rownames(score_df)))
    score_df
  })
  
  # Subset by common drugs
  common <- Reduce(intersect, lapply(score_list, rownames))
  score_list <- lapply(score_list, function(score_df) {
    score_df[rownames(score_df) %in% common, ]
  })
  
  # Add labels
  score_list <- lapply(score_list, function(score_df) {
    score_df$label <- ifelse(rownames(score_df) %in% approved_drugs,1,0)
    score_df
  })
  
  score_list
}

compute_ndcg <- function(
    scores  = NULL,
    disease = NULL
) {
  k_range <- c(5, 10, 15, nrow(scores$ReCAST))
  k_names <- k_range
  k_names[k_names == nrow(scores$ReCAST)] <- "Overall"
  
  ndcg_scores <- Map(function(score, m) {
    data.frame(
      score = sapply(k_range, function(k) {
        ndcg(score$label, score$drug_score, k)
      }),
      k = as.character(k_names),
      method = m
    )
  }, scores, names(scores))

  ndcg_scores <- do.call(rbind, ndcg_scores)
  ndcg_scores$method <- factor(ndcg_scores$method, levels = names(scores))
  ndcg_scores$k <- factor(ndcg_scores$k, levels = k_names)
  ndcg_scores$disease <- switch(disease,
                                "GBM" = "Glioblastoma",
                                "TNBC" = "Triple Negative Breast Cancer",
                                "MEL"  = "Melanoma"
  )
  ndcg_scores
}


plot_ndcg <- function(
    ndcg_scores = NULL,
    disease     = NULL
) {
  
  methods <- levels(ndcg_scores$method)
  colors <- c("#5AAE61", "#8172C6", "#C77C7C")

  labels <- sapply(seq_along(methods), function(i) {
    paste0("<span style='color:", colors[i], "'>", methods[i], "</span>")
  })
  
  title <- switch(disease,
                  "GBM" = "Glioblastoma",
                  "TNBC" = "Triple Negative Breast Cancer",
                  "MEL"  = "Melanoma"
  )
  
  ggplot(ndcg_scores, aes(x = k, y = score, fill = method)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9) +
    scale_fill_manual(values = colors, labels = labels, breaks = unique(ndcg_scores$method)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.25))) +
    theme_minimal() +
    labs(
      title = title,
      x = "Top k",
      y = "NDCG",
      fill = NULL
    ) +
    guides(fill = guide_legend(override.aes = list(fill = NA, color = NA))) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(
        size = 24,
        face = "bold",
        hjust = 0.5
      ),
      
      # Add a visible square
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_text(size = 22, color = "black"),
      axis.title = element_text(size = 22),
      
      # Legend styling
      legend.position = c(0.5, 0.9),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = ggtext::element_markdown(size = 20, face = "bold"),
      legend.key.height = unit(1.6, "lines")
    ) +
    geom_text(aes(label = round(score, 2)),
              position = position_dodge(0.9),
              vjust = -0.3,
              size = 6) 
}

##--------------------------------- GBM -------------------------------------
diseases <- c("TNBC", "GBM", "MEL")
results <- lapply(diseases, function(disease) {
  
  approved_drugs <- load_approved(disease)
  overall_scores <- load_scores(
    score_path = "data", 
    disease = disease, 
    approved_drugs = approved_drugs, 
    type = "overall"
  )
  
  ndcg_scores <- compute_ndcg(overall_scores, disease)
  plot <- plot_ndcg(ndcg_scores, disease)
  
  return(list(ndcg_scores = ndcg_scores,
              plot = plot
              )
         )
})
names(results) <- diseases

ndcg_scores <- do.call(rbind, lapply(results, "[[", "ndcg_scores"))
plots <- lapply(results, "[[", "plot")

wrap_plots(plots, ncol = 3)
ggsave(
  filename = file.path("data", "figures", "overall_ndcg_plot.png"), 
  plot = last_plot(), 
  width = 24, 
  height = 8, 
  units = "in", 
  dpi = 600
)
