#' @title Build CMap drug reference profiles
#'
#' @description
#' For each requested biological context (tissue or cell line), extracts the
#' relevant drug response signatures from the LINCS L1000 data
#' (GSE70138 and GSE92742), computes expression matrices, and writes three
#' files per context: '*_exprMatrix.txt', '*_gene_info.txt' and '*_drug_info.txt'.
#' If the output files for a context already exist, they are loaded directly
#' instead of being recomputed. After all contexts are processed, the function
#' merges them into a single reference object.
#'
#' @param sig_data_paths A named list of paths to the LINCS signature data, including
#' the cell metadata file, gene metadata file, GSE70138 signature metadata file,
#' GSE92742 signature metadata file, GSE70138 GCTX expression file and GSE92742
#' GCTX expression file.
#' @param contexts Character vector of biological context names to process
#' (e.g., c("MCF7", "HL60"), or c("breast")).
#' @param context_type Biological context type, either 'tissue' or 'cell'.
#' @param cores Number of cores to be used (1 by default).
#' @param out_path Output directory path (created if absent).
#'
#' @return A named list with: 'drug_info', drug_expr_matrix' and 'gene_info'.
#' Output files are written to 'out_path'.
#' @export

generate_reference <- function(sig_data_paths,
                               contexts,
                               context_type = "tissue",
                               cores        = 1,
                               out_path     = ".") {

  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # Read metadata
  cell_meta       <- read.table(sig_data_paths$cell_info_path, sep = "\t", header = TRUE, quote = "")
  gene_meta       <- read.table(sig_data_paths$gene_info_path, sep = "\t", header = TRUE, quote = "")
  sig_meta_70138  <- read.delim(sig_data_paths$gse70138_sig_info_path, sep = "\t", stringsAsFactors = FALSE)
  sig_meta_92742  <- read.delim(sig_data_paths$gse92742_sig_info_path, sep = "\t", stringsAsFactors = FALSE)

  # Select metadata column based on context
  if (context_type == "tissue") {
    valid_contexts <- unique(as.character(cell_meta$primary_site))
    valid_contexts <- valid_contexts[valid_contexts != "-666"]
    col_name   <- "primary_site"
  } else {
    valid_contexts <- cell_meta$cell_id[cell_meta$primary_site != "-666"]
    col_name   <- "cell_id"
  }

  selected_contexts <- intersect(contexts, valid_contexts)
  discarded_contexts <- setdiff(contexts, valid_contexts)
  if (is.null(contexts) || length(selected_contexts) == 0) {
    stop(paste0(
      paste(contexts, collapse = ", "),
      " not found. Available contexts: ",
      paste(valid_contexts, collapse = ", ")
    ))
  }

  if (length(discarded_contexts) > 0) {
    message(paste0(
        paste(discarded_contexts, collapse = ", "),
        " not found. Available contexts: ",
        paste(valid_contexts, collapse = ", ")
    ))
  }

  message(paste0(
    "Returning reference profiles for: ", paste(selected_contexts, collapse = ", ")
  ))

  # Returns valid sig_ids for a given set of cell lines, excluding replicates
  select_valid_sig_ids <- function(sig_meta, cell_lines, is_breast_context) {
    keep <- sig_meta$cell_id %in% cell_lines & sig_meta$pert_type == "trt_cp"
    if (is_breast_context) {
      keep <- keep & (sig_meta$pert_id != "BRD-K18910433")
    }
    sig_ids <- sig_meta$sig_id[keep]
    sig_ids[!grepl("REP\\.", sig_ids)]
  }

  # Parses a GCTX file for the given sig_ids and returns an expression matrix
  load_gctx <- function(gctx_path, sig_ids) {
    if (length(sig_ids) > 0) {
      expression_mat <- cmapR::parse_gctx(gctx_path, cid = sig_ids)@mat
      return(as.data.frame(expression_mat))
    } else {
      return(NULL)
    }
  }

  # Builds a drugs metadata dataframe from sig_meta
  build_drug_metadata <- function(sig_meta, standardized_sig_ids) {
    sig_meta$sig_id <- gsub(":", "_", sig_meta$sig_id)
    drug_info_df <- data.frame(
      instance_id       = sig_meta$sig_id,
      cmap_name         = paste(sig_meta$pert_iname, sig_meta$pert_id, sep = "_"),
      concentration     = sig_meta$pert_idose,
      duration          = sig_meta$pert_itime,
      cell              = sig_meta$cell_id,
      catalog_name      = sig_meta$pert_id,
      treatment         = paste0(sig_meta$pert_iname, "_", sig_meta$sig_id),
      stringsAsFactors  = FALSE
    )
    subset(drug_info_df, instance_id %in% standardized_sig_ids)
  }

  read_context <- function(ctx_name) {
    gene_info   <- read.table(paste0(out_path, "/", ctx_name, "_gene_info.txt"),
                              sep = "\t", header = TRUE, quote = "")
    drug_info   <- read.table(paste0(out_path, "/", ctx_name, "_drug_info.txt"),
                              sep = "\t", header = TRUE, quote = "")
    expr_matrix <- readRDS(paste0(out_path, "/", ctx_name, "_exprMatrix.rds"))
    list(drug_info = drug_info, expr_matrix = expr_matrix, gene_info = gene_info)
  }

  # Per context processing
  process_one_context <- function(context_name) {

    context_file <- gsub(" ", "-", context_name)
    expr_matrix_file <- paste0(context_file, "_exprMatrix.rds")
    gene_info_file   <- paste0(context_file, "_gene_info.txt")
    drug_info_file   <- paste0(context_file, "_drug_info.txt")

    if (file.exists(file.path(out_path, expr_matrix_file))) {
      return(read_context(context_file))
    }

    cell_lines  <- cell_meta$cell_id[cell_meta[[col_name]] == context_name]
    is_breast   <- (context_name == "breast")

    sig_ids_70138 <- select_valid_sig_ids(sig_meta_70138, cell_lines, is_breast)
    sig_ids_92742 <- select_valid_sig_ids(sig_meta_92742, cell_lines, is_breast)

    expr_mat_70138 <- load_gctx(sig_data_paths$gse70138_gctx_path, sig_ids_70138)
    expr_mat_92742 <- load_gctx(sig_data_paths$gse92742_gctx_path, sig_ids_92742)

    has_70138 <- !is.null(expr_mat_70138)
    has_92742 <- !is.null(expr_mat_92742)
    if (!has_70138 && !has_92742) {
      message(paste0(
        "No instances found for context: ",
        context_name
      ))
    }

    if (has_70138 && has_92742) {
      expr_mat_92742 <- expr_mat_92742[rownames(expr_mat_70138), , drop = FALSE]
      combined_expr_mat <- cbind(expr_mat_70138, expr_mat_92742)
    } else if (has_70138) {
      combined_expr_mat <- expr_mat_70138
    } else {
      combined_expr_mat <- expr_mat_92742
    }

    # Standardize column names
    standardized_sig_ids        <- gsub(":", "_", colnames(combined_expr_mat))
    colnames(combined_expr_mat) <- standardized_sig_ids

    # Replace string sig_ids with integer indexes
    integer_col_indices         <- as.character(seq_len(ncol(combined_expr_mat)))
    colnames(combined_expr_mat) <- integer_col_indices

    # -- expr matrix --
    saveRDS(combined_expr_mat,
            file  = file.path(out_path, expr_matrix_file))

    # -- gene info --
    gene_info_output <- gene_meta[, 1:2]
    colnames(gene_info_output) <- c("ID", "Gene.Symbol")
    write.table(gene_info_output,
                file  = file.path(out_path, gene_info_file),
                quote = FALSE, row.names = FALSE, sep = "\t")

    # -- drug info --
    drug_info_combined <- rbind(
      build_drug_metadata(sig_meta_70138, standardized_sig_ids),
      build_drug_metadata(sig_meta_92742, standardized_sig_ids)
    )
    drug_info_combined$instance_id  <- seq_len(nrow(drug_info_combined))
    write.table(drug_info_combined,
                file  = file.path(out_path, drug_info_file),
                quote = FALSE, row.names = FALSE, sep = "\t")

    return(list(drug_info = drug_info_combined,
                expr_matrix = combined_expr_mat,
                gene_info = gene_info_output))
  }

  # Entry point
  n_cores  <- min(cores, parallel::detectCores() - 1, length(selected_contexts))
  raw_data <- parallel::mclapply(selected_contexts, process_one_context,
                                 mc.cores = max(1, n_cores))

  names(raw_data) <- selected_contexts

  # Merge contexts sequentially
  offset <- 0
  drug_info_per_ctx  <- setNames(vector("list", length(selected_contexts)), selected_contexts)
  expr_mat_per_ctx   <- setNames(vector("list", length(selected_contexts)), selected_contexts)

  for (ctx_name in selected_contexts) {
    ctx_data    <- raw_data[[ctx_name]]
    gene_info   <- ctx_data$gene_info
    drug_info   <- ctx_data$drug_info
    expr_matrix <- ctx_data$expr_matrix

    # Shift instance IDs
    drug_info$instance_id <- drug_info$instance_id + offset
    colnames(expr_matrix)  <- as.character(as.integer(colnames(expr_matrix)) + offset)

    # Reorder rows and assign gene symbols
    expr_matrix <- expr_matrix[as.character(gene_info$ID), , drop = FALSE]
    rownames(expr_matrix) <- gene_info$Gene.Symbol

    drug_info_per_ctx[[ctx_name]] <- drug_info
    expr_mat_per_ctx[[ctx_name]]  <- expr_matrix
    offset <- offset + nrow(drug_info)
  }

  # Build result
  return(
    list(
      drug_info        = do.call(rbind, drug_info_per_ctx),
      drug_expr_matrix = do.call(cbind, unname(expr_mat_per_ctx)),
      gene_info        = gene_info
    )
  )
}
