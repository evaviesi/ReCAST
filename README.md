# ReCAST: Repurposing of drugs via a Cell cluster- and Action-aware Scoring Tool

ReCAST prioritizes drugs across cell populations and by individual cell type 
accounting for protein-protein interactions (PPIs) and drug mechanisms of action (MoAs).

## Installation
First, install devtools and the required packages:
```r
install.packages("devtools")
install.packages("Seurat")
install.packages("effsize")
install.packages("rlang")
install.packages("qs2")
install.packages("PRROC")
install.packages("R.utils")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(c("SingleR", "cmapR", "celldex"))

devtools::install_github('immunogenomics/presto')
devtools::install_github("lanagarmire/Asgard")
devtools::install_github("lhbcb/scDrugLink")
``` 
 
You can install ReCAST from GitHub using:
```r
devtools::install_github("InfOmics/ReCAST")
```

Loading:

```r
library("ReCAST")
```

## Usage

### **1. Load scRNA-seq data**
As an example, we use the melanoma (`MEL`) dataset to illustrate the steps 
required to reproduce the drug repurposing pipeline and the findings presented 
in the ReCAST article.

The `MEL_preprocess.R` script provides the code for data download, filtering and 
preprocessing, while the `ReCAST_score.R` script contains the code to compute 
the final drug score.

```r
library("Seurat")

out_path <- "out_data"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

tissue <- c("skin")
disease <- "MEL"

data(SC_obj)
```

#### **1.1. Single-cell processing, alignment and annotation** 
**Note**: The example dataset (`SC_obj`) has undergone quality control (QC), 
but still needs to be processed, aligned and annotated before proceeding to the 
next steps. If you already have a processed and annotated dataset, you can skip 
this section and proceed directly to the differential expression analysis.

```r
source("utils/align_cells.R")

SC_list <- SplitObject(SC_obj, split.by = "sample")
SC_obj <- align_cells(SC_list = SC_list,
                      CellCycle = TRUE,
                      anchor_features = 2000,
                      ndim = 15,
                      npcs = 15,
                      resolution = 0.4,
                      by_CellType = TRUE)
                             
DefaultAssay(SC_obj) <- "RNA"
SC_obj <- JoinLayers(SC_obj)

```

### **2. Perform differential gene expression analysis across cell types**
For each cell type, this step performs differential expression (DE) analysis using 
the Seurat Wilcoxon test to identify marker genes. The DE results are stored in a list.

```r
source("utils/filter_data.R")
source("utils/get_de_genes.R")

DefaultAssay(SC_obj) <- "RNA"

# Case and control sample names
case    <- unique(grep("MOF", SC_obj@meta.data$sample, value=TRUE))
control <- unique(grep("Normal", SC_obj@meta.data$sample, value=TRUE))

## Add "type" column to metadata
SC_obj@meta.data$type <- "Case"
SC_obj@meta.data[SC_obj@meta.data$sample %in% control,]$type <- "Control"

# Filter cell types with fewer than 10 cells per sample
annotation <- "celltype"
SC_obj <- filter_data(SC_obj     = SC_obj,
                      annotation = annotation,
                      min_cells  = 10)

# Get DE genes
deg_genes <- get_de_genes(SC_obj     = SC_obj,
                          annotation = annotation)

```

### **3. Load drug signatures from L1000 CMap**  
Drug perturbation signatures are downloaded from the Connectivity Map (CMap) L1000 
perturbational profiles GSE70138 and GSE92742 available in GEO. These signatures 
are then used to generate context-specific drug references with the 
`generate_reference` function. 

All drugs retrieved by `generate_reference` are stored and subsequently used as 
input for the DrugBank parsing step to retrieve their corresponding drug targets.

**Note**: For the `MEL` dataset, we selected `skin` as the biological context and 
`tissue` as the context type. 

```r
source("utils/download_signatures.R")

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
write.csv(all_drugs, file.path(out_path, "all_drugs_skin.csv"))
```

#### **3.1. Filter drugs having at least k instances** 
**Note**: We selected drugs with at least 2 instances to ensure reliable drug-level 
statistics and reduce sensitivity to single measurements, but the minimum number of 
instances can be customized.

```r
k <- 2
instance_count <- table(drug_ref_profiles$drug_info$cmap_name)
instance_count <- instance_count[instance_count >= k]
drug_ref_profiles$drug_info <- drug_ref_profiles$drug_info[
  drug_ref_profiles$drug_info$cmap_name %in% names(instance_count), 
]
drug_ref_profiles$drug_expr_matrix <- drug_ref_profiles$drug_expr_matrix[
  , as.character(drug_ref_profiles$drug_info$instance_id)
]
keep_drugs <- intersect(gsub("_.*","",drug_ref_profiles$drug_info$cmap_name), FDA_drugs) 
```

#### **3.2. Prepare drug target PPIs and drug-target MoAs**
To generate the `drug_targets_moa` object, containing the mechanisms of action
associated with drug targets downloaded from [DrugBank](https://go.drugbank.com/releases),
the data are processed using an XML parser and data extraction procedures described 
in the `DrugBank_parser.ipynb` notebook. The drugs listed in the `all_drugs_skin` 
file are used as input, with the `interactions` flag set to `True`.
 
The `drug_targets_ppi` data are generated from  protein-protein interactions of 
drug targets downloaded from [STRING](https://string-db.org/). These data are first 
generated by extracting drug-target interactions using the `DrugBank_parser.ipynb` 
notebook with the `interactions` flag set to `False`. The target sets are then
extended using STRING PPIs, as described in the `PPI_extended_targets.R` script.

```r
# Load PPI and MoA
drug_targets_moa <- read.csv(file.path(out_path, "drug_targets_moa_skin.csv"))
drug_targets_ppi <- read.csv(file.path(out_path, "drug_targets_ppi_950_skin.csv"))

# Filter PPI and MoA
drug_targets_ppi <- drug_targets_ppi[tolower(drug_targets_ppi$drug_name) %in% keep_drugs, ]
drug_targets_moa <- drug_targets_moa[tolower(drug_targets_moa$drug_name) %in% keep_drugs, ]
```

### **4. Compute drug statistics** 
For each cell type, this step computes drug statistics (p-value and FDR) using 
gene expression signature matching and the Kolmogorov-Smirnov test (adapted 
from the [`Asgard`](https://github.com/lanagarmire/Asgard) approach).

**Note**: Here, the `concordance` parameter is set to `negative` for 
disease-reversing drug perturbational profiles, but it can also be set to `positive` 
for disease-mimicking ones. The `drug_type` parameter is set to `FDA`, but can 
also be set to `compounds` or `ALL`. Important genes associated with drugs can 
optionally be returned.

```r
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
```

### **5. Compute drug action weight** 
This step builds the Drug2Cell (D2C) matrix using the R implementation designed 
by [`scDrugLink`](https://github.com/LHBCB/scDrugLink), following the Python 
implementation developed by the Teichmann Lab (https://github.com/Teichlab/drug2cell).

It then computes the drug action weight for each cell type, which combines the 
FDR and Cliff's delta from the Wilcoxon rank-sum test with an additional weight 
based on the drug's mechanism of action (inhibition or activation), determined 
according to the `concordance` parameter.

```r
# Get common gene list (CMap and SC data)
gene_list_common <- Reduce(intersect,
                           list("seurat" = rownames(SC_obj@assays$RNA),
                                "drug"   = drug_ref_profiles$gene_info$Gene.Symbol))

# Build drug2cell matrix
d2c_mat <- ReCAST::compute_d2c_matrix(SC_obj        = SC_obj,
                                      common_genes  = gene_list_common,
                                      target_ppi_df = drug_targets_ppi,
                                      out_path      = out_path)

# Compute drug action weight
drug_action_weight <- ReCAST::compute_drug_action(SC_obj                 = SC_obj,
                                                  d2c_mat                = d2c_mat,
                                                  cluster_degs           = deg_genes,
                                                  target_interactions_df = drug_targets_moa,
                                                  annotation             = annotation,
                                                  concordance            = concordance,
                                                  cores                  = 1,
                                                  out_path               = out_path)

```

### **6. Compute final drug score** 
The overall drug score is computed by summing the weighted contributions of each 
cell type, which combine drug statistics and drug action weights (as described in 
the ReCAST article), resulting in a ranking of repurposed drugs.

```r
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
```

#### **6.1. Evaluate results**
To evaluate the results in terms of AUC and AUPRC, the set of positive drugs for 
the corresponding disease is loaded from the `positive_drugs` object provided in 
the `ReCAST` package.

```r
source("utils/evaluate.R")

# Load positive drugs
positive_drugs <- rownames(positive_drugs[positive_drugs$MEL == 1, ])
evaluate(drug_score_df = drug_score, positive_drugs = positive_drugs)
```

## Citation


## Acknowledgement
ReCAST is derived from the [Asgard](https://github.com/lanagarmire/Asgard) and [scDrugLink](https://github.com/LHBCB/scDrugLink) packages.
We thank their authors for making these tools available. See the `LICENSE` file for the applicable license terms.

