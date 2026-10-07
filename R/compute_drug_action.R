#' @title Compute drug action weights per cluster/cell type
#'
#' @description
#' For each cell type, quantifies how strongly each drug shifts the distribution
#' of its target genes activity scores between diseased and healthy cells.
#'
#' The computation is parallelized across clusters using 'mclapply()'.
#'
#' @param SC_obj Seurat object with a metadata column matching 'annotation'.
#' @param d2c_mat Drug-2-cell matrix from 'compute_d2c_matrix' function.
#' @param cluster_degs Named list with a data.frame of DE genes per cluster.
#' @param target_interactions_df A data.frame with columns 'drug_name',
#' 'gene_names' and 'actions'.
#' @param annotation Name of the metadata column with cell type labels.
#' @param concordance Direction of concordance: either 'negative' or 'positive'.
#' @param cores Number of cores used to run the process in parallel.
#' @param seed Random seed for reproducibility (default 42).
#' @param out_path Output directory path.
#' @param out_name Output file name.
#'
#' @return A data.frame of drug weights (clusters/cell types x drugs).
#' @export

compute_drug_action <- function(SC_obj,
                                d2c_mat,
                                cluster_degs,
                                target_interactions_df,
                                annotation,
                                concordance  = "negative",
                                cores        = 1,
                                seed         = 42,
                                out_path     = "./",
                                out_name     = "action_weight") {

  set.seed(seed)

  Idents(SC_obj)  <- annotation
  cell_annot      <- as.character(unique(Idents(SC_obj)))

  # Process by cell type
  process_cell_annot <- function(cell_annot_name) {
    SC_obj_subset     <- subset(SC_obj, idents = cell_annot_name)
    disease_mask      <- as.integer(SC_obj_subset@meta.data$type == "Case")
    d2c_mat_subset    <- d2c_mat[, colnames(SC_obj_subset), drop = FALSE]

    diseased_cell_ids <- colnames(d2c_mat_subset)[disease_mask == 1]
    healthy_cell_ids  <- colnames(d2c_mat_subset)[disease_mask == 0]
    diseased_d2c      <- d2c_mat_subset[, diseased_cell_ids, drop = FALSE]
    healthy_d2c       <- d2c_mat_subset[, healthy_cell_ids,  drop = FALSE]

    cluster_deg  <- cluster_degs[[cell_annot_name]]
    drug_names   <- rownames(diseased_d2c)

    compute_drug_stats <- function(drug_name) {
      diseased_scores  <- as.numeric(diseased_d2c[drug_name, ])
      healthy_scores   <- as.numeric(healthy_d2c[drug_name,  ])

      if (var(diseased_scores) == 0 || var(healthy_scores) == 0) {
        warning(drug_name, ": zero variance in at least one group; skipping tests.")
        return(list(p_value = 1, cliff_delta = 0, action_weight = 0))
      }

      wilcox_pval     <- wilcox.test(diseased_scores, healthy_scores,
                                   alternative = "two.sided", paired = FALSE,
                                   conf.int = TRUE, conf.level = 0.95)$p.value
      cliff_delta_val <- effsize::cliff.delta(diseased_scores, healthy_scores)$estimate

      # MoA weight
      action_weight <- 0
      if (drug_name %in% target_interactions_df$drug_name) {
        drug_row       <- target_interactions_df[target_interactions_df$drug_name == drug_name, ]
        drug_genes     <- strsplit(drug_row$gene_names, ";")[[1]]
        drug_actions   <- strsplit(drug_row$actions,    ";")[[1]]
        targets_in_deg <- drug_genes %in% rownames(cluster_deg)

        if (sum(targets_in_deg) > 0) {
          drug_genes   <- drug_genes[targets_in_deg]
          drug_actions <- drug_actions[targets_in_deg]

          if (concordance == "negative") {
            action_sign   <- ifelse(drug_actions == "inhibitor",  1, -1)
            action_mask   <- action_sign == 1
          } else {
            action_sign   <- ifelse(drug_actions == "activator",  1, -1)
            action_mask   <- action_sign == 1
          }

          target_scores <- cluster_deg[drug_genes, "score"]
          moa_weight <- sum(as.integer(target_scores[action_mask]) * as.integer(action_sign[action_mask]) > 0)
          action_weight <- moa_weight / sum(targets_in_deg)
        }
      }

      list(p_value = wilcox_pval, cliff_delta = cliff_delta_val, action_weight = action_weight)
    }

    stats_by_drug <- lapply(drug_names, compute_drug_stats)
    names(stats_by_drug) <- drug_names

    p_values       <- sapply(stats_by_drug, `[[`, "p_value") + 1e-6
    cliff_deltas   <- sapply(stats_by_drug, `[[`, "cliff_delta")
    action_weights <- sapply(stats_by_drug, `[[`, "action_weight")

    adj_p_values <- p.adjust(p_values, method = "BH")

    # Final action weight
    weights <- exp(cliff_deltas * (-log10(adj_p_values))) * exp(action_weights)
    list(weights = weights)
  }

  # Parallel computation over cell types
  effective_cores <- max(1, min(cores, parallel::detectCores() - 1, length(cell_annot)))
  cell_annot_results <- parallel::mclapply(
    cell_annot, process_cell_annot,
    mc.cores = effective_cores
  )
  names(cell_annot_results) <- cell_annot

  weight_list    <- lapply(cell_annot_results, `[[`, "weights")

  drug_names <- tolower(rownames(d2c_mat))

  cell_annot_df <- function(cell_annot_vec) {
    df <- as.data.frame(t(cell_annot_vec))
    colnames(df) <- drug_names
    df
  }

  weight_df <- do.call(rbind, lapply(weight_list, cell_annot_df))
  rownames(weight_df) <- cell_annot

  write.csv(weight_df, file.path(out_path, paste0(out_name,  ".csv")), row.names = TRUE)

  return(weight_df)
}
