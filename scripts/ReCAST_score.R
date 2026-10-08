# Load packages and utilities --------------------------------------------------

library(ReCAST)
library(qs2)
library(Seurat)
source("utils/filter_data.R")
source("utils/get_de_genes.R")
source("utils/download_signatures.R")
source("utils/evaluate.R")


# Set parameters ---------------------------------------------------------------

set.seed(42)
out_path  <- "out_data"  # set output directory path
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

## Melanoma
tissue  <- c("skin")
disease <- "MEL"

## TNBC
# tissue  <- c("breast")
# disease <- "TNBC"

## GBM
# tissue  <- c("central nervous system")
# disease <- "GBM"


# Load data --------------------------------------------------------------------

## Load pre-processed scRNA-seq and set "RNA" as the default assay
SC_obj <- qs_read(file.path(out_path, "SC_obj.qs2"))
DefaultAssay(SC_obj) <- "RNA"
SC_obj <- JoinLayers(SC_obj)


# Differential expression analysis  --------------------------------------------

## Case and control sample names
case    <- unique(grep("MOF", SC_obj@meta.data$sample, value=TRUE))
control <- unique(grep("Normal", SC_obj@meta.data$sample, value=TRUE))

## Add "type" column to metadata
SC_obj@meta.data$type <- "Case"
SC_obj@meta.data[SC_obj@meta.data$sample %in% control,]$type <- "Control"

## Filter cell types with fewer than 10 cells per sample
annotation <- "celltype"
SC_obj <- filter_data(SC_obj     = SC_obj,
                      annotation = annotation,
                      min_cells  = 10)

## Get DE genes
deg_genes <- get_de_genes(SC_obj     = SC_obj,
                          annotation = annotation)


# Generate tissue-specific drug reference profiles -----------------------------

ref_path <- file.path(out_path, "L1000_CMap")
download_signatures(ref_path)

sig_paths <- list(cell_info_path         = file.path(ref_path, "GSE92742_Broad_LINCS_cell_info.txt"),
                  gene_info_path         = file.path(ref_path, "GSE92742_Broad_LINCS_gene_info.txt"),
                  gse70138_sig_info_path = file.path(ref_path, "GSE70138_Broad_LINCS_sig_info.txt"),
                  gse92742_sig_info_path = file.path(ref_path, "GSE92742_Broad_LINCS_sig_info.txt"),
                  gse70138_gctx_path     = file.path(ref_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"),
                  gse92742_gctx_path     = file.path(ref_path, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"))

drug_ref_profiles <- ReCAST::generate_reference(sig_data_paths = sig_paths,
                                                contexts       = tissue,
                                                context_type   = "tissue",
                                                cores          = 1,
                                                out_path       = file.path(out_path, "DrugReference"))

all_drugs <- unique(drug_ref_profiles$drug_info$cmap_name)
write.csv(all_drugs, file.path(out_path, paste0("all_drugs_", gsub(" ", "_", tissue), ".csv")))

# Use `all_drugs` as input to the `DrugBank_parser.ipynb` notebook to generate:
#   - `drug_targets_moa_{tissue}` with `interactions` set to `True`.
#   - `drug_targets_{tissue}` with `interactions` set to `False`. Then pass this file
#     to `PPI_extended_targets.R` to generate `drug_targets_ppi_950_{tissue}`.

## Load drug target PPIs and drug-target MoAs
drug_targets_ppi <- read.csv(paste0(out_path, "/drug_targets_ppi_950_", gsub(" ", "_", tissue), ".csv"))
drug_targets_moa <- read.csv(paste0(out_path, "/drug_targets_moa_", gsub(" ", "_", tissue), ".csv"))

## Select drugs with at least 2 instances
k <- 2
instance_count <- table(drug_ref_profiles$drug_info$cmap_name)
instance_count <- instance_count[instance_count >= k]
drug_ref_profiles$drug_info <- drug_ref_profiles$drug_info[drug_ref_profiles$drug_info$cmap_name %in% names(instance_count), ]
drug_ref_profiles$drug_expr_matrix <- drug_ref_profiles$drug_expr_matrix[, as.character(drug_ref_profiles$drug_info$instance_id)]
keep_drugs <- intersect(gsub("_.*","",drug_ref_profiles$drug_info$cmap_name), FDA_drugs)

# Compute drug statistics ------------------------------------------------------

drug_type <- "FDA"
imp_genes <- TRUE
concordance <- "negative"  # or "positive"

drug_stats <- ReCAST::compute_drug_statistics(cluster_degs      = deg_genes,
                                              drug_ref_profiles = drug_ref_profiles,
                                              repurposing_unit  = "drug",
                                              CEG_threshold     = 0.05,
                                              concordance       = concordance,
                                              drug_type         = drug_type,
                                              imp_genes         = imp_genes,
                                              cores             = 1)


# Compute drug action weight ---------------------------------------------------

## Filter drug target PPIs and drug-target MoAs
drug_targets_ppi <- drug_targets_ppi[tolower(drug_targets_ppi$drug_name) %in% keep_drugs, ]
drug_targets_moa <- drug_targets_moa[tolower(drug_targets_moa$drug_name) %in% keep_drugs, ]

## Get common gene list (scRNA-seq and CMap)
gene_list_common <- Reduce(intersect,
                           list("seurat" = rownames(SC_obj@assays$RNA),
                                "drug"   = drug_ref_profiles$gene_info$Gene.Symbol))

## Build drug2cell matrix
d2c_mat <- ReCAST::compute_d2c_matrix(SC_obj        = SC_obj,
                                      common_genes  = gene_list_common,
                                      target_ppi_df = drug_targets_ppi,
                                      out_path      = out_path)

drug_action_weight <- ReCAST::compute_drug_action(SC_obj                 = SC_obj,
                                                  d2c_mat                = d2c_mat,
                                                  cluster_degs           = deg_genes,
                                                  target_interactions_df = drug_targets_moa,
                                                  annotation             = annotation,
                                                  concordance            = concordance,
                                                  cores                  = 1,
                                                  out_path               = out_path)


# Compute final drug score -----------------------------------------------------

cell_metadata <- SC_obj@meta.data

drug_score <- ReCAST::compute_drug_score(cell_metadata  = cell_metadata,
                                         annotation     = annotation,
                                         case           = case,
                                         cluster_degs   = deg_genes,
                                         sig_data_paths = sig_paths,
                                         score_drugs    = drug_stats,
                                         action_weights = drug_action_weight,
                                         contexts       = tissue,
                                         context_type   = "tissue",
                                         clusters       = NULL,
                                         concordance    = concordance,
                                         drug_type      = drug_type,
                                         imp_genes      = imp_genes,
                                         imp_genes_thr  = 0.1,
                                         cores          = 1,
                                         out_path       = out_path)


# Evaluate results -------------------------------------------------------------

##  Load positive drugs
positive_drugs <- rownames(positive_drugs[positive_drugs$MEL == 1, ])
evaluate(drug_score_df = drug_score, positive_drugs = positive_drugs)

