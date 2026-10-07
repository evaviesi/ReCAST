#' @title Compute drug scores across cell clusters
#'
#' @description
#' Integrates drug statistics, cluster proportions, drug weights, and the
#' CMap gene expression profiles to produce a single score for each drug.
#'
#' The computation is parallelized across drugs using 'mclapply()'.
#'
#' @param cell_metadata A data.frame of cell metadata with column 'sample'
#' and a column matching 'annotation'.
#' @param annotation Name of the metadata column with cell type labels.
#' @param case Case sample names to be evaluated.
#' @param cluster_degs Named list with a data.frame of differentially
#' expressed genes per cluster.
#' @param sig_data_paths A named list of paths to the LINCS signature data, including
#' the cell metadata file, gene metadata file, GSE70138 signature metadata file,
#' GSE92742 signature metadata file, GSE70138 GCTX expression file and GSE92742
#' GCTX expression file.
#' @param score_drugs Output of 'compute_drug_statistics' function.
#' @param action_weights Output of 'compute_drug_action' function.
#' @param clusters Subset of cluster names to be evaluated.
#' @param contexts Character vector of biological context names to process
#' (e.g., c("MCF7", "HL60"), or c("breast")).
#' @param context_type Biological context type, either 'tissue' or 'cell'.
#' @param concordance Direction of concordance: either 'negative' or 'positive'.
#' @param drug_type One of 'FDA', 'compounds', or 'ALL'.
#' @param imp_genes Logical; whether to return important genes.
#' @param imp_genes_thr Percentage of top important genes to return.
#' @param cores Number of cores to be used (1 by default).
#' @param out_path Output directory path.
#'
#' @return A data.frame with columns named 'drug_score', 'p_val', 'fdr'.
#' When imp_genes = TRUE: additional columns with gene names for each cluster.
#' @export

compute_drug_score <- function(cell_metadata,
                               annotation,
                               case,
                               cluster_degs,
                               sig_data_paths,
                               score_drugs,
                               action_weights,
                               contexts,
                               context_type  = "tissue",
                               clusters      = NULL,
                               concordance   = "negative",
                               drug_type     = "FDA",
                               imp_genes     = FALSE,
                               imp_genes_thr = 0.1,
                               cores         = 1,
                               out_path      = "./") {

  # Input validation
  if (!(drug_type %in% c("FDA", "compounds", "ALL"))) {
    stop("drug_type must be 'FDA', 'compounds', or 'ALL'.")
  }

  if (imp_genes_thr <= 0 || imp_genes_thr > 1) {
    stop("imp_genes_thr must be in range ]0,1]")
  }

  # Load important genes (if available)
  if (imp_genes) {
    cluster_drug_stats  <- score_drugs$drug
    cluster_gene_pvals  <- score_drugs$pvalue
  } else {
    cluster_drug_stats  <- score_drugs
    cluster_gene_pvals  <- NULL
  }

  # Subset selected clusters
  if (length(clusters) > 0) {
    clusters <- intersect(as.character(clusters), unique(cell_metadata[[annotation]]))
    if (length(clusters) == 0) {
      stop("No valid clusters found. Available: ",
           paste(unique(cell_metadata[[annotation]]), collapse = ", "))
    }
    cell_metadata       <- subset(cell_metadata, cell_metadata[[annotation]] %in% clusters)
    cluster_drug_stats  <- cluster_drug_stats[clusters]
    cluster_degs        <- cluster_degs[clusters]
  }

  # Cluster proportions in case samples
  if (length(case) > 0) {
    cell_metadata <- subset(cell_metadata, sample %in% case)
  } else {
    stop("No case samples found.")
  }

  cell_annot            <- cell_metadata[[annotation]]
  cluster_size_table    <- table(cell_annot)
  cluster_proportions   <- round(100 * cluster_size_table / nrow(cell_metadata), 2)

  drug_cluster_df <- data.frame()

  for (cluster_name in names(cluster_drug_stats)) {
    drug_stats  <- cluster_drug_stats[[cluster_name]]

    # Filter reference profiles
    candidate_drugs <- drug_stats$drug_name
    if (drug_type != "ALL") {
      if (drug_type == "FDA") {
        candidate_drugs <- intersect(drug_stats$drug_name, FDA_drugs)
      } else {
        candidate_drugs <- setdiff(drug_stats$drug_name, FDA_drugs)
      }
    }
     if (length(candidate_drugs) > 0) {
      drug_stats      <- subset(drug_stats, drug_name %in% candidate_drugs)
      drug_cluster_df <- rbind(drug_cluster_df, data.frame(
        drug_name    = drug_stats$drug_name,
        drug_id      = drug_stats$drug_id,
        cluster      = cluster_name,
        cluster_prop = cluster_proportions[cluster_name],
        p_value      = drug_stats$p_val,
        fdr          = drug_stats$fdr,
        row.names    = NULL,
        stringsAsFactors = FALSE
      ))
    }
  }
  drug_cluster_df <- unique(drug_cluster_df)
  drug_cluster_df$weighted_prop <- drug_cluster_df$cluster_prop * (-log10(drug_cluster_df$fdr))

  if (anyNA(drug_cluster_df)) {
    View(drug_list)
    stop("NA values detected in drug_cluster_df.")
  }

  all_drug_names <- unique(drug_cluster_df$drug_id)

  # Combined p-value across clusters for each drug
  if (length(unique(names(cluster_drug_stats))) > 1) {
    combined_pvalues <- tapply(drug_cluster_df$p_value, drug_cluster_df$drug_id, compute_chisq_pvalues)
  } else {
    combined_pvalues <- setNames(drug_cluster_df$p_value, drug_cluster_df$drug_id)
  }

  # Load CMap expression profiles ----------------------------------------
  # Read metadata
  cell_meta       <- read.table(sig_data_paths$cell_info_path, sep = "\t", header = TRUE, quote = "")
  gene_meta       <- read.table(sig_data_paths$gene_info_path, sep = "\t", header = TRUE, quote = "")
  sig_meta_70138  <- read.delim(sig_data_paths$gse70138_sig_info_path, sep = "\t", stringsAsFactors = FALSE)
  sig_meta_92742  <- read.delim(sig_data_paths$gse92742_sig_info_path, sep = "\t", stringsAsFactors = FALSE)

  if (context_type == "cell") {
    cell_lines <- subset(cell_meta, cell_id %in% contexts)$cell_id
  } else {
    cell_lines <- subset(cell_meta, primary_site %in% contexts)$cell_id
  }

  load_gctx <- function(col_meta, gctx_path) {
    sig_idx  <- which(col_meta$cell_id %in% cell_lines & col_meta$pert_id %in% all_drug_names)
    sig_ids  <- col_meta$sig_id[sig_idx]
    if (length(sig_ids) == 0) return(list(exprs = NULL, meta = NULL))

    meta <- col_meta[, c("sig_id", "pert_id")]
    rownames(meta) <- meta$sig_id
    meta <- meta[sig_ids, ]

    exprs_mat  <- as.data.frame(cmapR::parse_gctx(gctx_path, cid = sig_ids)@mat)
    exprs_mat$gene_id <- rownames(exprs_mat)
    merged_mat <- merge(exprs_mat, gene_meta, by.x = "gene_id", by.y = "pr_gene_id")
    exprs <- merged_mat[, c("pr_gene_symbol", sig_ids)]

    list(exprs = exprs, meta = meta, sig_ids = sig_ids)
  }

  data_92742 <- load_gctx(sig_meta_92742, sig_data_paths$gse92742_gctx_path)
  data_70138 <- load_gctx(sig_meta_70138, sig_data_paths$gse70138_gctx_path)

  has_92742 <- !is.null(data_92742$exprs)
  has_70138 <- !is.null(data_70138$exprs)

  if (!has_70138 && !has_92742) stop("No expression data found for the requested contexts.")

  if (has_70138 && has_92742) {
    drug_exprs <- merge(data_92742$exprs, data_70138$exprs, by = "pr_gene_symbol")
    rownames(drug_exprs) <- drug_exprs[, "pr_gene_symbol"]
    drug_exprs[["pr_gene_symbol"]] <- NULL
    drug_metadata  <- rbind(data_92742$meta, data_70138$meta)
  } else if (has_70138) {
    drug_exprs <- data_70138$exprs
    rownames(drug_exprs) <- drug_exprs[, "pr_gene_symbol"]
    drug_exprs[["pr_gene_symbol"]] <- NULL
    drug_metadata  <- data_70138$meta
  } else {
    drug_exprs <- data_92742$exprs
    rownames(drug_exprs) <- drug_exprs[, "pr_gene_symbol"]
    drug_exprs[["pr_gene_symbol"]] <- NULL
    drug_metadata  <- data_92742$meta
  }

  colnames(action_weights) <- tolower(gsub("\\.", "-", colnames(action_weights)))

  # ============================================================
  # Parallel drug scoring
  # ============================================================
  message("Computing drug scores ...")

  score_by_drug <- function(drug_name) {
    treatment_sig_ids    <- subset(drug_metadata, pert_id == drug_name)$sig_id
    drug_mean_expression <- rowMeans(drug_exprs[, treatment_sig_ids, drop = FALSE])

    drug_scores  <- drug_cluster_df[drug_cluster_df$drug_id == drug_name, ]
    total_score <- 0
    cluster_score <- setNames(numeric(length(cluster_degs)), names(cluster_degs))
    drug_lower <- unique(tolower(drug_scores[["drug_name"]]))

    if (drug_lower %in% colnames(action_weights)) {
      drug_weight <- action_weights[, drug_lower, drop = FALSE]
    } else {
      drug_weight <- data.frame(matrix(1, nrow = nrow(action_weights), ncol = 1,
                              dimnames = list(rownames(action_weights), drug_lower)))
    }

    for (cluster_name in names(cluster_degs)) {
      cluster_row <- drug_scores[drug_scores$cluster == cluster_name, ]

      sig_degs <- cluster_degs[[cluster_name]]
      sig_degs <- sig_degs[sig_degs$adj_p_val <= 0.05, , drop = FALSE]

      treatable_genes <- intersect(rownames(sig_degs), names(drug_mean_expression))
      if (length(treatable_genes) == 0) next

      deg_score        <- sig_degs[treatable_genes, "score"]
      drug_concordance <- -deg_score * drug_mean_expression[treatable_genes]

      if (concordance == "negative") {
        treated_genes <- drug_concordance[drug_concordance > 0]
      } else if (concordance == "positive") {
        treated_genes <- drug_concordance[drug_concordance < 0]
      } else {
        stop("concordance must be 'negative' or 'positive'.")
      }

      treated_fraction  <- length(treated_genes) / length(treatable_genes)
      contribution      <- (cluster_row$cluster_prop / 100) * (-log10(cluster_row$fdr)) * treated_fraction * drug_weight[cluster_name, ]

      total_score  <- total_score + contribution
      cluster_score[cluster_name] <- contribution
    }

    list(score = total_score, cluster_contributions = cluster_score)
  }

  effective_cores  <- max(1, min(cores, parallel::detectCores()- 1, length(all_drug_names)))
  res_by_drug <- parallel::mclapply(all_drug_names, score_by_drug, mc.cores = effective_cores)
  names(res_by_drug) <- all_drug_names

  drug_scores_res          <- sapply(res_by_drug, `[[`, "score")
  cluster_contribution_mat <- as.data.frame(t(sapply(res_by_drug, `[[`, "cluster_contributions")))

  scored_drugs_df <- data.frame(
    drug_score = drug_scores_res,
    p_val      = combined_pvalues[all_drug_names],
    fdr        = p.adjust(combined_pvalues[all_drug_names], method = "BH"),
    row.names  = all_drug_names
  )
  scored_drugs_df <- scored_drugs_df[order(scored_drugs_df$drug_score, decreasing = TRUE), ]

  # Compute important genes by drug and cluster (optional)
  if (imp_genes) {
    if (is.null(cluster_gene_pvals)) stop("imp_genes requires pvalue!")
    message("Computing important genes ...")

    cluster_filtered_pvals <- lapply(names(cluster_drug_stats), function(cluster_name) {
      pmin_df <- cluster_gene_pvals[["aggregated_pmin"]][[cluster_name]]
      pmax_df <- cluster_gene_pvals[["aggregated_pmax"]][[cluster_name]]
      min_gene_pval <- pmin(pmin_df, pmax_df)

      # Filtering for DE genes p-value
      sig_gene_names <- rownames(cluster_degs[[cluster_name]])[
        cluster_degs[[cluster_name]]$adj_p_val <= 0.05
      ]
      min_gene_pval[rownames(min_gene_pval) %in% sig_gene_names, , drop = FALSE]
    })
    names(cluster_filtered_pvals) <- names(cluster_drug_stats)

    for (cluster_name in names(cluster_filtered_pvals)) {
      cluster_pval <- cluster_filtered_pvals[[cluster_name]]
      important_genes <- sapply(rownames(scored_drugs_df), function(drug_name) {
        # Filtering for pmin/pmax p-value
        drug_gene_pvals <- cluster_pval[, drug_name, drop = FALSE]
        significant     <- subset(drug_gene_pvals, drug_gene_pvals <= 0.05)
        top_genes       <- significant[significant <= quantile(significant, imp_genes_thr, na.rm = TRUE), , drop = FALSE]
        paste(rownames(top_genes), collapse = ";")

      })

      scored_drugs_df[[paste0("Cluster_", cluster_name, "_genes")]] <- important_genes
    }
  }

  # Remove duplicated drugs
  drug_map <- setNames(drug_stats$drug_name, drug_stats$drug_id)
  scored_drugs_df$drug_name <- drug_map[rownames(scored_drugs_df)]
  scored_drugs_df <- scored_drugs_df[!duplicated(scored_drugs_df$drug_name), ]
  scored_drugs_df$drug_id <- rownames(scored_drugs_df)
  rownames(scored_drugs_df) <- scored_drugs_df$drug_name
  scored_drugs_df$drug_name <- NULL

  cluster_contribution_mat <- cluster_contribution_mat[rownames(cluster_contribution_mat) %in% scored_drugs_df$drug_id, ]
  rownames(cluster_contribution_mat) <- drug_map[rownames(cluster_contribution_mat)]

  # Saving results
  file_prefix <- if (concordance == "negative") "/neg_" else "/pos_"
  write.csv(scored_drugs_df,          paste0(out_path, file_prefix, "drug_scores.csv"),    quote = FALSE)
  write.csv(cluster_contribution_mat, paste0(out_path, file_prefix, "cluster_scores.csv"), quote = FALSE)

  return(scored_drugs_df)
}

