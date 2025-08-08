#!/usr/bin/env Rscript

library(biomaRt)

# Construct archived Ensembl BioMart host URL
# For example, version = "104" ==> host = "https://apr2021.archive.ensembl.org"
# Ensembl archived URLs correspond to release date, approximate mapping required:
# See https://www.ensembl.org/info/website/archives/index.html

# Some example archive dates mapping to versions:
version_map <- list(
  "84" = "mar2016",
  "85" = "jul2016",
  "86" = "oct2016",
  "87" = "dec2016",
  "88" = "mar2017",
  "89" = "may2017",
  "90" = "aug2017",
  "91" = "dec2017",
  "92" = "apr2018",
  "93" = "jul2018",
  "94" = "oct2018",
  "95" = "jan2019",
  "96" = "apr2019",
  "97" = "jul2019",
  "98" = "sep2019",
  "99" = "jan2020",
  "100" = "apr2020",
  "101" = "aug2020",
  "102" = "nov2020",
  "103" = "feb2021",
  "104" = "may2021",
  "105" = "dec2021",
  "106" = "apr2022",
  "107" = "jul2022",
  "108" = "oct2022",
  "109" = "feb2023",
  "110" = "jul2023",
  "111" = "jan2024",
  "112" = "may2024",
  "113" = "oct2024",
  "114" = "may2025"
)


fetch_unique_genes <- function(ensembl_version, output_file) {
  
  if (!(ensembl_version %in% names(version_map))) {
    stop("Version not supported or invalid. Please add the mapping for your Ensembl version in version_map.")
  }
  
  host_url <- paste0("https://", version_map[[as.character(ensembl_version)]], ".archive.ensembl.org")

  ensembl <- useEnsembl(biomart = "genes",
                       dataset = "hsapiens_gene_ensembl")


  # Connect to the archive BioMart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl",
                  host=host_url)

  # Query attributes: ensembl_gene_id, external_gene_name, gene_biotype
  genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                 mart = mart)

  # Remove duplicates by ensembl_gene_id
  genes_unique <- genes[!duplicated(genes$ensembl_gene_id), ]

  # Write to output file
  write.table(genes_unique, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = c("Ensembl_ID", "Gene_Name", "Gene_Biotype"))

  message(sprintf("Unique genes for Ensembl version %s saved to %s", ensembl_version, output_file))
}

# Example of running the function:
# fetch_unique_genes("104", "genes_ens104.tsv")

# To enable running from command line:
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    cat("Usage: Rscript fetch_ensembl_genes.R <ensembl_version> <output_file>\n")
    quit(status = 1)
  }
  version <- args[1]
  outfile <- args[2]
  fetch_unique_genes(version, outfile)
}
