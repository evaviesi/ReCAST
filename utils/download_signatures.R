library(R.utils)

# ------------------------------------------------------------------------------
# This function downloads the Connectivity Map (CMap) signatures from the GEO
# portal, extracts and renames the compressed files.
# ------------------------------------------------------------------------------

download_signatures <- function(data_path) {

  if (!dir.exists(data_path)) {
    dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Download Phase I LINCS L1000 CMap data (GEO: GSE92742)
  system2("wget", args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2EMODZ%5Fn473647x12328%2Egctx%2Egz"))
  system2("wget", args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz"))
  system2("wget", args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fcell%5Finfo%2Etxt%2Egz"))
  system2("wget", args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%2Etxt%2Egz"))

  # Download Phase II LINCS L1000 CMap data (GEO: GSE70138)
  system2("wget", args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%5Fn118050x12328%5F2017%2D03%2D06%2Egctx%2Egz"))
  system2("wget",  args = c("-c", "-P", data_path, "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5Fsig%5Finfo%5F2017%2D03%2D06%2Etxt%2Egz"))

  # Unzip files
  gz_files <- list.files(data_path, full.names = TRUE)
  for (gz_file in gz_files) {
    gunzip(gz_file, remove = TRUE)
  }

  # Rename files
  file.rename(
    file.path(data_path, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"),
    file.path(data_path, "GSE70138_Broad_LINCS_sig_info.txt")
  )

  file.rename(
    file.path(data_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"),
    file.path(data_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
  )

}

