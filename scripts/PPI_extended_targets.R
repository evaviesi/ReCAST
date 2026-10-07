# Libraries ---------------------------------------------------------------
library(data.table)
library(gprofiler2)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Input -------------------------------------------------------------------
data_path <- "data"
ppi_path <- file.path("PPI", "9606.protein.links.detailed.v12.0.txt.gz")
exp_score = 950
ppi_url <- "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz"
disease <- "MEL"
tissue <- "skin"
drug_targets_path <- file.path("data", disease, "input", "targets_skin.csv")

# Main --------------------------------------------------------------------

# Download PPI from string-db
dir.create("PPI")
options(timeout = max(3600, getOption("timeout")))
download.file(
  url = ppi_url,
  destfile = ppi_path
)

# Load and filter PPI
ppi_full <- fread(ppi_path, data.table = F)
ppi <- ppi_full[ppi_full$experimental >= exp_score, ]

# Remove 9606 prefix from proteins ENSP
ppi$protein1 <- gsub("9606.", "", ppi$protein1)
ppi$protein2 <- gsub("9606.", "", ppi$protein2)

# Convert protein ENSP to gene ENSG
ppi$protein1 <- gconvert(
  ppi$protein1,
  organism   = "hsapiens",
  target     = "ENSG",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na  = FALSE
)
ppi$protein2 <- gconvert(
  ppi$protein2,
  organism   = "hsapiens",
  target     = "ENSG",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na  = FALSE
)

# Keep only ENSG and experimental columns, removing rows with any na
protein1 <- ppi$protein1$target
protein2 <- ppi$protein2$target
score <- ppi$experimental
ppi_mapped <- data.frame(protein1 = protein1, protein2 = protein2, combined_score = score)
ppi_mapped <- ppi_mapped[!(apply((is.na(ppi_mapped)),1,any)), ]

# Convert gene ENSG to gene Symbol, removing rows with any na
ppi_mapped$protein1  <- mapIds(
  org.Hs.eg.db, 
  keys      = ppi_mapped$protein1, 
  column    = 'SYMBOL', 
  keytype   = 'ENSEMBL', 
  multiVals = "first"
)
ppi_mapped$protein2  <- mapIds(
  org.Hs.eg.db, 
  keys      = ppi_mapped$protein2, 
  column    = 'SYMBOL', 
  keytype   = 'ENSEMBL', 
  multiVals = "first")
ppi_mapped <- ppi_mapped[!(apply((is.na(ppi_mapped)),1,any)), ]


# Extend drug targets with 1-hop target neighborhood
drug_targets <- read.csv(drug_targets_path)

drug_targets$gene_names <- vapply(seq_len(nrow(drug_targets)), function(i) {
  # Extract targets of current drug
  ppi_genes <- unlist(strsplit(drug_targets$gene_names[i], ";"))
  
  # Collect 1-hop neighborhood PPI genes
  near_genes <- c(
    ppi_mapped$protein2[ppi_mapped$protein1 %in% ppi_genes],
    ppi_mapped$protein1[ppi_mapped$protein2 %in% ppi_genes]
  )
  
  # Combine 1-hop neighborhood genes with direct targets
  all_genes <- unique(c(ppi_genes, near_genes))
  paste0(all_genes, collapse = ";")
}, character(1))

# Save extended drug-targets interactions
write.csv(
  drug_targets, 
  file.path(data_path, disease, "output", 
            paste0("targets_ppi_", exp_score, "_", tissue, ".csv")),
  row.names = FALSE
)
