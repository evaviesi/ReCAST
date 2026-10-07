# Libraries and Packages --------------------------------------------------
library(ReCAST)
library(ROCR)
library(PRROC)
library(ggplot2)
library(patchwork)
# Function ----------------------------------------------------------------

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
      ASGARD     = file.path("ASGARD",     paste0(disease, "_drug_scores.csv")),
      scDrugLink = file.path("scDrugLink", paste0(disease, "_drug_scores.csv")),
      ReCAST     = file.path("output", "neg_drug_scores.csv")
    )
  } else {
    score_files <- list(
      ASGARD     = file.path("ASGARD",     paste0(disease, "_individual_cell_type_drug_scores.csv")),
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

# Plot either ROC or PR curve
metric_plot <- function(
    metric_df = NULL,
    metric    = "PR",
    title     = NULL,
    labels    = NULL,
    methods   = NULL,
    area      = NULL,
    rnd_value = NULL
) {
  
  # Add random line
  rnd_line <- if (metric == "ROC") {
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray"
    )
  } else {
    geom_hline(
      yintercept = rnd_value, 
      linetype = "dashed", 
      color = "gray"
    ) 
  }
  
  # Define input values
  if (metric == "ROC") {
    plot <- ggplot(metric_df, aes(x = FPR, y = TPR, color = method))
    area_type <- "AUC"
  } else if (metric == "PR") {
    plot <- ggplot(metric_df, aes(x = Recall, y = Precision, color = method))
    area_type <- "AUPRC"
  }
  
  # Add colors and labels
  colors <- c("#5AAE61", "#8172C6", "#C77C7C")
  labels <- sapply(seq_along(methods), function(i) {
    paste0(
      "<span style='color:", colors[i], "'>", methods[i], 
      " (", area_type, "= ", round(area[[methods[i]]], 2), ")</span>"
    )
  })
  
  plot <- plot +
    geom_line(linewidth = 1, alpha = 0.9) +
    rnd_line +
    scale_color_manual(values = colors, labels = labels, breaks = unique(metric_df$method)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = title,
      color = NULL) +
    guides(color = guide_legend(override.aes = list(linetype = 0, size = 0))) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      
      # Add a visible square
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18),
      
      # Legend styling
      legend.position = if (metric == "ROC") {c(0.67, 0.12)} else {c(0.65, 0.75)},
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = ggtext::element_markdown(size = 16, face = "bold")
    )
      
  plot
}


# Compute and save ROC curves plot with AUROC values. 
# Returns the ROC dataframe.
ROC_AUC <- function(
    scores       = NULL,
    figure_dir   = "data/figures",
    disease      = NULL,
    by_cluster   = NULL,
    cluster      = NULL,
    cluster_prop = NULL
) {
  
  methods <- names(scores)
  
  # Defne and create figures output dir
  figure_path <- file.path(figure_dir, disease)
  dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
  
  # Make prediction
  pred <- lapply(scores, function(score_df) {
    prediction(score_df$drug_score, score_df$label)
  })
  
  # ROC - AUC
  roc <- lapply(pred, function(p) {
    performance(p, "tpr", "fpr")
  })
  auc <- lapply(pred, function(p) {
    performance(p, "auc")@y.values[[1]]
  })

  # Create a dataframe of ROC curves
  roc_df <- do.call(rbind, Map(function(obj, m) {
    data.frame(FPR = obj@x.values[[1]], TPR = obj@y.values[[1]], method = m)
  }, roc, names(roc)))
  
  # Add cell type and disease column for excel
  if (by_cluster) {
    roc_df$celltype <- cluster
  }
  roc_df$disease <- switch(disease,
    "GBM"  = "Glioblastoma",
    "TNBC" = "Triple Negative Breast Cancer",
    "MEL"  = "Melanoma"
  )
  
  # ROC-AUC plot
  if (!by_cluster) {
    title <- switch(disease, 
                   "GBM"  = "Glioblastoma",
                   "TNBC" = "Triple Negative Breast Cancer",
                   "MEL"  = "Melanoma"
    )
  } else {
    title <- paste0(cluster, " (", cluster_prop, "%)")
  }

  plot <- metric_plot(roc_df, "ROC", title, labels, methods, auc)
  
  # Plot and save results: write if overall score, append on plots otherwise
  if (!by_cluster) {
    print(plot)
    ggsave(
      filename = file.path(figure_path, "roc_plot.png"),
      plot = plot,
      width = 6,
      height = 6,
      units = "in",
      dpi = 600
    )
    return(roc_df)
  } else {
    return(list(metric_df = roc_df, plot = plot))
  }
}


# Compute and save PR curves plot with AUPRC values
# Returns the PR dataframe.
PR_AUPRC <- function(
  scores       = NULL,
  approved     = NULL,
  figure_dir   = "data/figures",
  disease      = NULL,
  by_cluster   = NULL,
  cluster      = NULL,
  cluster_prop = NULL
) {
  
  methods <- names(scores)
  
  # Define and create figures output dir
  figure_path <- file.path(figure_dir, disease)
  dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
  
  # Retrieve positive and negative scores
  pos <- lapply(scores, function(score_df) {
    score_df$drug_score[score_df$label == 1]
  })
  neg <- lapply(scores, function(score_df) {
    score_df$drug_score[score_df$label == 0]
  })
  
  # Compute PR and AUPRC
  pr <- Map(function(p,n) {
    pr.curve(p, n, curve = TRUE)
  }, pos, neg)
  auprc <- setNames(sapply(pr, "[", "auc.integral"), methods)
  
  # Create a dataframe of PR curves
  pr_df <- do.call(rbind, Map(function(obj, m) {
    data.frame(Recall = obj$curve[,1], Precision = obj$curve[,2], method = m)
  }, pr, names(pr)))
  
  # Add cell type and disease column for excel
  if (by_cluster) {
    pr_df$celltype <- cluster
  }
  pr_df$disease <- switch(disease,
    "GBM"  = "Glioblastoma",
    "TNBC" = "Triple Negative Breast Cancer",
    "MEL"  = "Melanoma"
  )
  
  # PR-AUPRC plot
  if (!by_cluster) {
    title <- switch(disease, 
                    "GBM"  = "Glioblastoma",
                    "TNBC" = "Triple Negative Breast Cancer",
                    "MEL"  = "Melanoma"
    )
  } else {
    title <- paste0(cluster, " (", cluster_prop, "%)")
  }
  
  # Add random line
  pos <- length(intersect(approved, rownames(scores$ReCAST)))
  neg <- nrow(scores$ReCAST) - pos
  rnd_line <- pos / (pos + neg)
  
  plot <- metric_plot(pr_df, "PR", title, labels, methods, auprc, rnd_line)
  
  # Plot and save results: write if overall score, append on plots otherwise
  if (!by_cluster) {
    print(plot)
    ggsave(
      filename = file.path(figure_path, "pr_plot.png"), 
      plot = plot, 
      width = 6, 
      height = 6, 
      units = "in", 
      dpi = 600
    )
    return(pr_df)
  } else {
    return(list(metric_df = pr_df, plot = plot))
  }
}

# Pipeline for cluster-level ROC/PR curves. 
# Returns the ROC/PR dataframe containing all cluster metrics.
cluster_plot <- function(
    cluster_scores = NULL,
    cluster_prop   = NULL,
    disease        = NULL,
    figure_dir     = "data/figures",
    metric         = "PR",
    approved_drugs = NULL
) {
  
  # Define and create figures output dir
  figure_path <- file.path(figure_dir, disease)
  dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
  
  # Organize scores by cell line, keeping separate data frames for each tool
  pivot_scores <- lapply(names(cluster_scores[[1]])[-ncol(cluster_scores[[1]])], function(cell_line) {
    setNames(
      lapply(cluster_scores, function(df) {
        data.frame(
          drug_score = df[[cell_line]],
          label = df$label,
          row.names = rownames(df)
        )
      }),
      names(cluster_scores)
    )
  })
  names(pivot_scores) <- names(cluster_scores[[1]])[-ncol(cluster_scores[[1]])]
  
  # Computer ROC or PR for each cluster
  results <- lapply(names(pivot_scores), function(cell_line) {
    if (metric == "ROC") {
      ROC_AUC(
        scores       = pivot_scores[[cell_line]], 
        disease      = disease,
        by_cluster   = TRUE,
        cluster      = cell_line,
        cluster_prop = cluster_prop[[cell_line]]
      )
    } else {
      PR_AUPRC(
        scores       = pivot_scores[[cell_line]],
        disease      = disease,
        approved     = approved_drugs,
        by_cluster   = TRUE,
        cluster      = cell_line,
        cluster_prop = cluster_prop[[cell_line]]
      )
    }

  })
  names(results) <- names(pivot_scores)
  
  # Reorganize into common ROC dataframe and plot list
  metric_df <- do.call(rbind, lapply(results, "[[", "metric_df"))
  plots <- lapply(results, "[[", "plot")
  
  # Define plot parameters
  title = switch(disease, 
                 "GBM"  = "Glioblastoma",
                 "TNBC" = "Triple Negative Breast Cancer",
                 "MEL"  = "Melanoma"
  )
  n_col = switch(disease, 
    "GBM"  = 3,
    2
  )
  
  # Wraps cluster plots into unique plot and save it
  filename = ifelse(metric == "ROC", "roc_celltype_plot.png", "pr_celltype_plot.png")
  w_plot <- wrap_plots(plots, ncol = n_col) +  plot_annotation(title = title) &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  print(w_plot)
  ggsave(
    filename = file.path(figure_path, filename), 
    plot = w_plot, 
    width = 18, 
    height = 18, 
    units = "in", 
    dpi = 600
  )
  
  metric_df
}


# GBM ---------------------------------------------------------------------

disease = "GBM"
approved_drugs <- load_approved(disease)

# Overall score ROC and PR
gbm_overall_scores <- load_scores(
  score_path = "data", 
  disease = disease, 
  approved_drugs = approved_drugs, 
  type = "overall"
)
gbm_roc_df <- ROC_AUC(
  scores = gbm_overall_scores, 
  figure_dir = "data/figures", 
  disease = disease,
  by_cluster = FALSE
)
gbm_pr_df <- PR_AUPRC(
  scores = gbm_overall_scores, 
  figure_dir = "data/figures", 
  approved = approved_drugs,
  disease,
  by_cluster = FALSE
)

# Cluster-level score ROC and PR
gbm_cluster_scores <- load_scores(
  score_path = "data", 
  disease, 
  approved_drugs, 
  "cluster"
)
cluster_prop   <- readRDS(file.path("data", disease, "output", "cluster_proportions.rds"))
colnames(cluster_prop) <- trimws(colnames(cluster_prop))
cluster_prop <- cluster_prop[, intersect(colnames(gbm_cluster_scores$ReCAST), colnames(cluster_prop))]

gbm_cluster_roc_df <- cluster_plot(
  cluster_scores = gbm_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "ROC"
)

gbm_cluster_pr_df <- cluster_plot(
  cluster_scores = gbm_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "PR",
  approved_drugs = approved_drugs
)



# TNBC --------------------------------------------------------------------

disease = "TNBC"
approved_drugs <- load_approved(disease, "data")

# Overall score ROC and PR
tnbc_overall_scores <- load_scores(
  score_path = "data", 
  disease = disease, 
  approved_drugs = approved_drugs, 
  type = "overall"
)
tnbc_roc_df <- ROC_AUC(
  scores = tnbc_overall_scores, 
  figure_dir = "data/figures", 
  disease = disease,
  by_cluster = FALSE
)
tnbc_pr_df <- PR_AUPRC(
  scores = tnbc_overall_scores, 
  figure_dir = "data/figures", 
  approved = approved_drugs,
  disease,
  by_cluster = FALSE
)

# Cluster-level score ROC and PR
tnbc_cluster_scores <- load_scores(
  score_path = "data", 
  disease, 
  approved_drugs, 
  "cluster"
)
cluster_prop   <- readRDS(file.path("data", disease, "output", "cluster_proportions.rds"))
colnames(cluster_prop) <- trimws(colnames(cluster_prop))
cluster_prop <- cluster_prop[, intersect(colnames(tnbc_cluster_scores$ReCAST), colnames(cluster_prop))]

tnbc_cluster_roc_df <- cluster_plot(
  cluster_scores = tnbc_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "ROC"
)

tnbc_cluster_pr_df <- cluster_plot(
  cluster_scores = tnbc_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "PR",
  approved_drugs = approved_drugs
)



# MEL ---------------------------------------------------------------------
disease = "MEL"
approved_drugs <- load_approved(disease, "data")

# Overall score ROC and PR
mel_overall_scores <- load_scores(
  score_path = "data", 
  disease = disease, 
  approved_drugs = approved_drugs, 
  type = "overall"
)
mel_roc_df <- ROC_AUC(
  scores = mel_overall_scores, 
  figure_dir = "data/figures", 
  disease = disease,
  by_cluster = FALSE
)
mel_pr_df <- PR_AUPRC(
  scores = mel_overall_scores, 
  figure_dir = "data/figures", 
  approved = approved_drugs,
  disease,
  by_cluster = FALSE
)

# Cluster-level score ROC and PR
mel_cluster_scores <- load_scores(
  score_path = "data", 
  disease, 
  approved_drugs, 
  "cluster"
)
cluster_prop   <- readRDS(file.path("data", disease, "output", "cluster_proportions.rds"))
colnames(cluster_prop) <- trimws(colnames(cluster_prop))
cluster_prop <- cluster_prop[, intersect(colnames(mel_cluster_scores$ReCAST), colnames(cluster_prop))]

mel_cluster_roc_df <- cluster_plot(
  cluster_scores = mel_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "ROC"
)

mel_cluster_pr_df <- cluster_plot(
  cluster_scores = mel_cluster_scores, 
  cluster_prop   = cluster_prop, 
  disease        = disease,
  metric         = "PR",
  approved_drugs = approved_drugs
)
