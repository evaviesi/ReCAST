# Load packages and utilities --------------------------------------------------

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

## TNBC
# tissue  <- "breast"

## GBM
# tissue  <- "central-nervous-system"


# Generate tissue-specific drug reference profiles -----------------------------

ref_path <- file.path(out_path, "L1000_CMap")
download_signatures(ref_path)
dir.create(file.path(out_path, "AsgardDrugReference/"), recursive = TRUE, showWarnings = FALSE)

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
                         Output.Dir        = file.path(out_path, "AsgardDrugReference/"))


## Select drugs with at least 2 instances
drug_info <- read.table(file = paste0(out_path, "/AsgardDrugReference/", tissue, "_drug_info.txt"), sep="\t", header = T, quote = "")
drug_rank_matrix <- read.table(paste0(out_path, "/AsgardDrugReference/", tissue, "_rankMatrix.txt"), row.names = 1, header = T, check.names = FALSE)

k <- 2
instance_count <- table(drug_info$cmap_name)
instance_count <- instance_count[instance_count >= k]
drug_info <- drug_info[drug_info$cmap_name %in% names(instance_count), ]
drug_rank_matrix <- drug_rank_matrix[, as.character(drug_info$instance_id)]
drug_info$instance_id <- 1:nrow(drug_info)
colnames(drug_rank_matrix) <- 1:nrow(drug_info)
keep_drugs <- intersect(gsub("_.*","",drug_info$cmap_name), FDA_drugs)

## Save results
write.table(drug_info, paste0(out_path, "/AsgardDrugReference/", tissue, "_drug_info_sub.txt"), sep="\t", quote = FALSE)
write.table(drug_rank_matrix, paste0(out_path, "/AsgardDrugReference/", tissue, "_rankMatrix_sub.txt"), sep="\t", quote = FALSE)


# Differential expression analysis  --------------------------------------------

## Load pre-processed scRNA-seq and set "RNA" as the default assay
SC_obj <- qread(file.path(out_path, "SC_obj.qs"))
DefaultAssay(SC_obj) <- "RNA"
SC_obj <- JoinLayers(SC_obj)

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


# Load tissue-specific drug reference profiles ---------------------------------

gene_info <- read.table(file=paste0(out_path, "/AsgardDrugReference/", tissue, "_gene_info.txt"), sep = "\t", header = T, quote = "")
drug_info <- read.table(file=paste0(out_path, "/AsgardDrugReference/", tissue, "_drug_info_sub.txt"), sep = "\t", header = T, quote = "")

drug.ref.profiles <- Asgard::GetDrugRef(drug.response.path = paste0(out_path, "/AsgardDrugReference/", tissue, "_rankMatrix_sub.txt"),
                                        probe.to.genes = gene_info,
                                        drug.info = drug_info)


# Mono-drug repurposing for every cell type ------------------------------------

deg_genes <- lapply(deg_genes, function(df) {
  colnames(df) <- c("score", "adj.P.Val", "P.Value"); return(df) })

Drug.ident.res <- Asgard::GetDrug(gene.data         = deg_genes,
                                  drug.ref.profiles = drug.ref.profiles,
                                  repurposing.unit  = "drug",
                                  CEG.threshold     = 0.05,
                                  connectivity      = "negative",
                                  drug.type         = "FDA")


# Calculate final drug score ---------------------------------------------------

gse92742_gctx_path <- file.path(ref_path, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
gse70138_gctx_path <- file.path(ref_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")

cell_metadata         <- SC_obj@meta.data
cell_metadata$cluster <- SC_obj@meta.data[[annotation]]

Drug.score <- Asgard::DrugScore(cell_metadata      = cell_metadata,
                                cluster_degs       = deg_genes,
                                cluster_drugs      = Drug.ident.res,
                                tissue             = tissue,
                                case               = case,
                                gse92742_gctx_path = gse92742_gctx_path,
                                gse70138_gctx_path = gse70138_gctx_path)


# Evaluate results -------------------------------------------------------------

##  Load positive drugs
positive_drugs <- rownames(positive_drugs[positive_drugs$MEL == 1, ])
colnames(Drug.score) <- c("drug_score", "p_val", "fdr")
evaluate(drug_score_df = Drug.score, positive_drugs = positive_drugs)

