# Load packages and utilities --------------------------------------------------

library(scDrugLink)
library(Asgard)
library(qs)
source("utils/filter_data.R")
source("utils/get_de_genes.R")
source("utils/download_signatures.R")
source("utils/evaluate.R")


# Set parameters ---------------------------------------------------------------

set.seed(42)
out_path  <- "out_data"  # set output directory path
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

## Melanoma
tissue  <- "skin"
disease <- "MEL"

## TNBC
# tissue  <- "breast"
# disease <- "TNBC"

## GBM
# tissue  <- "central-nervous-system"
# disease <- "GBM"

# Load data --------------------------------------------------------------------

## Load pre-processed scRNA-seq and set "RNA" as the default assay
SC_obj <- qread(file.path(out_path, "SC_obj.qs"))
DefaultAssay(SC_obj) <- "RNA"
SC_obj <- JoinLayers(SC_obj)


# Generate tissue-specific drug reference profiles -----------------------------

ref_path <- file.path(out_path, "L1000_CMap")
download_signatures(ref_path)
dir.create(file.path(out_path, "scDrugLinkDrugReference/"), recursive = TRUE, showWarnings = FALSE)

cell.info         <- file.path(file.path(ref_path, "GSE92742_Broad_LINCS_cell_info.txt"))
gene.info         <- file.path(file.path(ref_path, "GSE92742_Broad_LINCS_gene_info.txt"))
GSE70138.sig.info <- file.path(file.path(ref_path, "GSE70138_Broad_LINCS_sig_info.txt"))
GSE92742.sig.info <- file.path(file.path(ref_path, "GSE92742_Broad_LINCS_sig_info.txt"))
GSE70138.gctx     <- file.path(file.path(ref_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"))
GSE92742.gctx     <- file.path(file.path(ref_path, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"))

Asgard::PrepareReference(cell.info         = cell.info,
                         gene.info         = gene.info,
                         GSE70138.sig.info = GSE70138.sig.info,
                         GSE92742.sig.info = GSE92742.sig.info,
                         GSE70138.gctx     = GSE70138.gctx,
                         GSE92742.gctx     = GSE92742.gctx,
                         Output.Dir        = file.path(out_path, "scDrugLinkDrugReference/"))


## Select drugs with at least 2 instances
drug_info <- read.table(file = paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_drug_info.txt"), sep="\t", header = T, quote = "")
drug_rank_matrix <- read.table(paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_rankMatrix.txt"), row.names = 1, header = T, check.names = FALSE)

k <- 2
instance_count <- table(drug_info$cmap_name)
instance_count <- instance_count[instance_count >= k]
drug_info <- drug_info[drug_info$cmap_name %in% names(instance_count), ]
drug_rank_matrix <- drug_rank_matrix[, as.character(drug_info$instance_id)]
drug_info$instance_id <- 1:nrow(drug_info)
colnames(drug_rank_matrix) <- 1:nrow(drug_info)
keep_drugs <- intersect(gsub("_.*","",drug_info$cmap_name), FDA_drugs)

## Save results
write.table(drug_info, paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_drug_info_sub.txt"), sep="\t", quote = FALSE)
write.table(drug_rank_matrix, paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_rankMatrix_sub.txt"), sep="\t", quote = FALSE)


# Building Drug2Cell matrix ----------------------------------------------------

gene_info <- read.table(file=paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_gene_info.txt"), sep="\t", header = T, quote = "")
common_gene_list <- Reduce(intersect, list("seurat" = rownames(SC_obj@assays["RNA"]$RNA),
                                           "drug"   = gene_info$Gene.Symbol))

## Load drug targets
# Use `unique(drug_info$cmap_name)` as input to the `DrugBank_parser.ipynb` notebook to
# generate `drug_targets_{tissue}` setting `interactions` to `False`.
drug_targets <- read.csv(paste0(out_path, "/drug_targets_", gsub(" ", "_", tissue), ".csv"))
drug_targets$drug_name_lower <- tolower(drug_targets$drug_name)
drug_targets  <- drug_targets[drug_targets$drug_name_lower %in% keep_drugs, ]
drug_label_df <- drug_targets # this parameters is globally captured to not change build_drug_target_d2c function

## Convert Assay5 to Assay class to not change build_drug_target_d2c function
SC_obj[["RNA"]] <- tryCatch(
  expr = {
    as(object = SC_obj[["RNA"]], Class = "Assay")
  },
  error = function(e) {
    data_mat   <- GetAssayData(SC_obj[["RNA"]], layer = "data")
    rna_v4 <- CreateAssayObject(
      data   = data_mat
    )
  }
)

d2c_mat <- scDrugLink::build_drug_target_d2c(dat            = SC_obj,
                                             gene_list      = common_gene_list,
                                             drug_target_df = drug_targets,
                                             out_path       = out_path)


# Estimate drug promotion/inhibition effects -----------------------------------

## Add "disease" and "cell_type" columns to metadata
SC_obj$disease   <- ifelse(SC_obj$sample %in% case, disease, "healthy")
SC_obj$cell_type <- SC_obj[[annotation]]

drug_prom_inh_weight <- scDrugLink::compute_drug_prom_inh(dat      = SC_obj,
                                                          d2c_mat  = d2c_mat,
                                                          disease  = disease,
                                                          out_type = "weight",
                                                          out_path = out_path)


# Differential expression analysis  --------------------------------------------

## Case and control sample names
case    <- unique(grep("MOF", SC_obj@meta.data$sample, value=TRUE))
control <- unique(grep("Normal", SC_obj@meta.data$sample, value=TRUE))

## Add "type" column to metadata
SC_obj@meta.data$type <- "Case"
SC_obj@meta.data[SC_obj@meta.data$sample %in% control,]$type <- "Control"

## Filter cell types with fewer than 10 cells per sample
annotation <- "celltype"
SC_obj     <- filter_data(SC_obj = SC_obj,
                          annotation = annotation,
                          min_cells  = 10)

## Get DE genes
method    <- "Seurat"
deg_genes <- get_de_genes(SC_obj     = SC_obj,
                          annotation = annotation,
                          control    = control,
                          case       = case,
                          method     = method)



# Estimate drug sensitivity/resistance effects ---------------------------------

deg_genes <- lapply(deg_genes, function(df) {
  colnames(df) <- c("score", "adj.P.Val", "P.Value"); return(df) })

drug_info <- read.table(file=paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_drug_info_sub.txt"), sep="\t", header = T, quote = "")
perturbation_matrix_path <- paste0(out_path, "/scDrugLinkDrugReference/", gsub(" ", "-", tissue), "_rankMatrix_sub.txt")

drug_sens_res_pval <- scDrugLink::compute_drug_sens_res(disease                  = disease,
                                                        perturbation_matrix_path = perturbation_matrix_path,
                                                        gene_info                = gene_info,
                                                        drug_info                = drug_info,
                                                        deg_list                 = deg_genes,
                                                        out_path                 = out_path)

# Calculate final drug score ---------------------------------------------------

gse92742_gctx_path <- file.path(ref_path, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
gse70138_gctx_path <- file.path(ref_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")

drug_score <- scDrugLink::compute_scdruglink_score(dat                  = SC_obj,
                                                   deg_list             = deg_genes,
                                                   drug_sens_res_pval   = drug_sens_res_pval,
                                                   tissue               = tissue,
                                                   case                 = disease,
                                                   gse92742_gctx_path   = gse92742_gctx_path,
                                                   gse70138_gctx_path   = gse70138_gctx_path,
                                                   drug_prom_inh_weight = drug_prom_inh_weight,
                                                   out_path             = out_path)


# Evaluate results -------------------------------------------------------------

##  Load positive drugs
positive_drugs <- rownames(positive_drugs[positive_drugs$MEL == 1, ])
evaluate(drug_score_df = drug_score, positive_drugs = positive_drugs)
