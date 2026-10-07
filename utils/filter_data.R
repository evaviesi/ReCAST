# ------------------------------------------------------------------------------
# This function filters data based on the number of cells per sample, minimum
# cell proportion and cell ratio.
# ------------------------------------------------------------------------------

filter_data <- function(SC_obj,
                        annotation = "celltype",
                        min_cells  = 3,
                        min_cell_proportion = NULL,
                        min_cell_ratio = NULL) {

  n_cell_by_type <- table(SC_obj$sample, SC_obj@meta.data[[annotation]])
  keep_cells     <- colnames(n_cell_by_type)[apply(n_cell_by_type, 2, function(x) all(x >= min_cells))]
  SC_obj         <- subset(SC_obj, !!rlang::sym(annotation) %in% keep_cells)

  if (!is.null(min_cell_proportion)) {
    n_cells    <- table(SC_obj[[annotation]])
    tot_cells  <- sum(n_cells)
    prop_cells <- n_cells/tot_cells
    keep_cells <- names(prop_cells)[(prop_cells >= min_cell_proportion)]
    SC_obj     <- subset(SC_obj, !!rlang::sym(annotation) %in% keep_cells)
  }

  if (!is.null(min_cell_ratio)) {
    n_cell_by_type  <- table(SC_obj$type, SC_obj@meta.data[[annotation]])
    case_to_control <- n_cell_by_type["Case", ] / n_cell_by_type["Control", ]
    keep_cells      <- names(case_to_control)[(min_cell_ratio <= case_to_control & case_to_control <= 1/min_cell_ratio)]
    SC_obj          <- subset(SC_obj, !!rlang::sym(annotation) %in% keep_cells)
  }
  message("Remaining cell types: ", paste0(unique(SC_obj@meta.data[[annotation]]), sep = " " ))
  return(SC_obj)
}
