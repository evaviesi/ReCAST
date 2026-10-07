# Libraries and Packages --------------------------------------------------
library(ReCAST)
library(circlize)
library(ComplexHeatmap)
# Function ----------------------------------------------------------------

# Filter positive drugs for disease of interest
load_approved <- function(
    disease   = NULL
) {
  rownames(ReCAST::positive_drugs)[ReCAST::positive_drugs[[disease]] == 1]
}

# Filter score data by common perturbagens, and add labels
load_and_filter_scores <- function(
    score_path     = "data",
    disease        = NULL,
    n_show         = 15
) {
  # ---------- Define score folder and file names ----------
  data_folder <- file.path(score_path, disease, "output")
  overall_names <- c("neg_drug_scores.csv", "pos_drug_scores.csv") 
  cluster_names <- c("neg_cluster_scores.csv", "pos_cluster_scores.csv")
  
  
  # ---------- Load and filter overall data ----------
  overall_scores <- lapply(overall_names, function(p) {
    # Load overall score
    score <- read.csv(file.path(data_folder, p), row.names = 1)
    
    # Filter by significant drugs and top n_show scoring drugs
    score <- score[score$fdr < 0.05,]
    score <- score[order(score$drug_score, decreasing = TRUE),]
    score <- score[1:n_show, , drop = FALSE]
    
    # Cast to matrix and retain only "drug_score"
    as.matrix(score[,"drug_score", drop = FALSE])
  })
  names(overall_scores) <- gsub("_.*", "", overall_names)
  
  
  # ---------- Load and filter cluster data ----------
  cluster_scores <- lapply(cluster_names, function(p) {
    # Concordance to subset coherently with overall score
    concordance <- gsub("_.*", "", p)
    
    # Load cluster score
    score <- read.csv(file.path(data_folder, p), row.names = 1)
    
    # Subset w.r.t. overall score
    as.matrix(score[rownames(overall_scores[[concordance]]), , drop = FALSE])
  })
  names(cluster_scores) <- gsub("_.*", "", cluster_names)
  
  list(cluster_scores = cluster_scores, overall_scores = overall_scores)
}

# Load cluster proportions
load_cluster_prop <- function(
    path    = "data",
    scores  = NULL, 
    disease = NULL 
) {
  cluster_prop <- readRDS(file.path(path, disease, "output", "cluster_proportions.rds"))
  colnames(cluster_prop) <- trimws(colnames(cluster_prop))
  cluster_prop[, colnames(scores$cluster_scores$neg)]
}


# Create cluster legend
create_cluster_legend <- function(
  cluster_prop = NULL
) {
  
  # Define color palette mapped to cluster proportions
  cluster_colors <- colorRamp2(
    c(min(cluster_prop), 1, 2, max(cluster_prop)),
    c("#f2f2f2", "#d9d9d9", "#a6a6a6", "#595959")
  )
  
  # Compute custom legend breaks (min, mean, max) and map corresponding colors
  # (taking from the cluster_colors palette the colors corresponding to the breaks)
  legend_breaks <- as.numeric(
    formatC(c(min(cluster_prop), mean(cluster_prop), max(cluster_prop)), 
            format = "e", digits = 1)
  )
  legend_colors <- cluster_colors(legend_breaks)
  
  # Configure legend appearance, layout, and typography
  Legend(
    title          = "Cluster\nproportions",
    title_position = "topcenter",
    col_fun        = cluster_colors,
    direction      = "horizontal",
    at             = legend_breaks,
    break_dist     = 1,
    border         = "darkgray",
    legend_width   = unit(30, "mm"),
    title_gp       = gpar(fontsize = 12),
    labels_gp      = gpar(fontsize = 10)
  )
}


rank_heatmap <- function(
  scores         = NULL,
  approved_drugs = NULL,
  cluster_prop   = NULL,
  concordance    = "neg"
) {
  
  # Define general palette for heatmaps
  palettes <- list(
    # Negative concordance
    "neg" = list(
      "overall_scores" = c("#d0e4ff", "#5a88dc", "#1d4fb4"),
      "cluster_scores" = c("#e6f2ff", "#89bce6", "#2574b3")
    ),
    # Positive concordance
    "pos" = list(
      "overall_scores" = c("#ffdcd2", "#ec5e3a", "#8f2d1e"),
      "cluster_scores" = c("#fff0eb", "#f59a84", "#b04a33")
    ),
    # Cluster
    "cluster" = c("#f2f2f2", "#d9d9d9", "#a6a6a6", "#595959")
  )
  
  # Extract concordance scores
  scores <- lapply(scores, "[[", concordance)
  
  # Define bold labels for positive perturbagens
  row_font <- ifelse(
    rownames(scores$overall_scores) %in% approved_drugs,  
    "bold", "plain"
  )
  
  # Create left row annotation for the heatmap
  cluster_colors <- colorRamp2(
    c(min(cluster_prop), 1, 2, max(cluster_prop)),
    palettes[["cluster"]]
  )
  left_annot <- rowAnnotation(
    "Cluster Proportions" = anno_simple(
      x   = cluster_prop, 
      col = cluster_colors, 
      gp  = gpar(col = "gray", lwd = 1)
    ),
    annotation_name_gp   = gpar(fontsize = 12),
    show_annotation_name = FALSE,
    width                = unit(3.8, "mm")
  )
  
  # Create heatmaps
  heatmaps <- lapply(names(scores), function(n_score) {
    # Data
    score   <- scores[[n_score]]

    # Colors
    palette <- palettes[[concordance]][[n_score]]
    colors  <- colorRamp2(
      as.numeric(quantile(unique(as.vector(score)), seq(0, 1, by = 0.1))),
      c(colorRampPalette(c(palette[1], palette[2]))(10), palette[3])
    )
    
    # Labels and titles
    row_labels <- if (n_score == "overall_scores") {
      ""
    } else {
      gsub("_", " ", colnames(score))
    }

    title <- paste0(
      ifelse(n_score == "overall_scores", "Overall ", "Cluster "),
      ifelse(concordance == "neg", "inhibition-\n", "activation-\n"),
      "aware score"
    )

    # Heatmap parameters 
    ht_params <- list(
      matrix    = t(score),
      col       = colors,
      rect_gp   = gpar(col = "gray", lwd = 1),
      row_order      = colnames(score),
      row_labels     = row_labels,
      row_names_gp   = gpar(fontsize = 12),
      row_names_side = "left",
      column_order   = rownames(score),
      column_names_rot      = 45,
      column_names_centered = FALSE,
      column_names_gp       = gpar(fontsize = 12, fontface = row_font),
      heatmap_legend_param  = list(
        title          = title,
        title_position = "topcenter",
        col_fun        = colors,
        border         = "darkgray",
        at = as.numeric(formatC(c(min(score,  na.rm = TRUE),
                                  mean(score, na.rm = TRUE),
                                  max(score,  na.rm = TRUE)),
                                format = "e", digits = 1)),
        break_dist   = 1,
        direction    = "horizontal",
        legend_width = unit(30, "mm"),
        title_gp     = gpar(fontsize = 12),
        labels_gp    = gpar(fontsize = 10)
      ),
      cluster_columns = FALSE,
      cluster_rows    = FALSE,
      width  = unit(nrow(score) * 6, "mm"),
      height = unit(ncol(score) * 6, "mm")
    )
    
    # Heatmap parameters for cluster_scores only
    if (n_score == "cluster_scores") {
      ht_params$column_title    = " "
      ht_params$left_annotation = left_annot
      ht_params$row_title       = "Cell type"
      ht_params$row_title_side  = "left"
      ht_params$row_title_gp    = gpar(fontsize = 14)
      ht_params$column_labels   = rownames(score)
    }

    # Create Heatmap
    do.call(Heatmap, ht_params)
  })
  names(heatmaps) <- names(scores)
  heatmaps
}

create_legend <- function(
  heatmaps = NULL   
) {
  # Extract legends from heatmaps
  legends <- setNames(lapply(heatmaps, function(h) 
    {do.call(Legend, h@matrix_legend_param)}
  ), names(heatmaps))
  
  # Create empty legend for padding
  empty_legend <- Legend(
    labels = "",
    pch = NA,
    grid_width = unit(10, "mm"),
  )
  
  # Create Heatmap legend
  packLegend(
    empty_legend,
    packLegend(
      empty_legend,
      legend_clusters,
      empty_legend,
      legends$cluster_scores,
      empty_legend,
      legends$overall_scores,
      direction = "horizontal"
    ),
    direction = "vertical"
  )
}


# Main ------------------------------------------------------------------

disease     <- "TNBC"
concordance <- "neg"
approved_drugs <- load_approved(disease)
scores <- load_and_filter_scores("data", disease)
cluster_prop <- load_cluster_prop("data", scores, disease)

legend_clusters <- create_cluster_legend(cluster_prop)

heatmaps <- rank_heatmap(scores, approved_drugs, cluster_prop, concordance)
legends <- create_legend(heatmaps)
ht_list <- Reduce(`%v%`, heatmaps)

plot <- draw(
  object            = ht_list,
  column_title      = "Drug",
  column_title_side = "bottom",
  column_title_gp   = gpar(fontsize = 14),
  show_heatmap_legend    = FALSE,
  annotation_legend_list = list(legends),
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  padding = unit(c(10, 10, 10, 10), "mm"),
)
png(
  file.path("data", "figures", disease, paste0(concordance, "_heatmap.png")),
  width = 2300, 
  height = 1800, 
  res = 300
)
print(plot)
dev.off()
