#!/usr/bin/env Rscript

library(biomaRt)

Sys.setenv(SSL_CERT_FILE="/etc/ssl/certs/ca-certificates.crt")
Sys.setenv(CURL_CA_BUNDLE="/etc/ssl/certs/ca-certificates.crt")

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


fetch_unique_genes <- function(ensembl_version, output_file, merging_strat="first") {
  
  if (!(ensembl_version %in% names(version_map))) {
    stop("Version not supported or invalid. Please add the mapping for your Ensembl version in the version_map of this script. Current versions avaliable betwen 84-114")
  }
  
  host_url <- paste0("https://", version_map[[as.character(ensembl_version)]], ".archive.ensembl.org")

  # Connect to the archive BioMart
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl",
                  host=host_url)

  # Query attributes: ensembl_gene_id, external_gene_name, gene_biotype
  attr <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "source")
  genes <- getBM(attributes = attr,
                 mart = mart)

  # Remove duplicates by ensembl_gene_id
  genes_unique <- genes[!duplicated(genes$ensembl_gene_id), ]

  # Make sure the names are unique
  if (merging_strat == "first") {
    # Sort so the proper chromosomes and protein coding genes are first, this avoids wierd contigs getting the _
    
    # Define the chromosome order
    chrom_order <- c(as.character(1:22), "X", "Y", "MT")

    # Create a factor for chromosome_name with custom levels
    genes_unique$chromosome_name <- factor(
      genes_unique$chromosome_name, 
      levels = c(chrom_order, sort(setdiff(unique(genes_unique$chromosome_name), chrom_order)))
    )

    # Sort on biotype first
    biotype_order <- ifelse(genes_unique$gene_biotype == "protein_coding", 0, 1)

    # Sort the dataframe
    genes_unique <- genes_unique[order(genes_unique$chromosome_name, biotype_order), ]    
  
    genes_unique[,"final_gene_name"] <- genes_unique[,"hgnc_symbol"]
    genes_unique[is.na(genes_unique$final_gene_name), "final_gene_name"] <- "NO_NAME"
    genes_unique[genes_unique$final_gene_name == "", "final_gene_name"] <- "NO_NAME"
    
    # Dedup NO_NAME
    genes_unique[genes_unique$final_gene_name == "NO_NAME", "final_gene_name"] <- make.unique(genes_unique[genes_unique$final_gene_name == "NO_NAME", "final_gene_name"], sep="_")
    
    # Dedup any other gene names
    genes_unique[,"final_gene_name"] <- make.unique(genes_unique[,"final_gene_name"], sep="_")    
    
    # Sort on back on chromosome position
    genes_unique <- genes_unique[order(genes_unique$chromosome_name, genes_unique$start_position), ]    
    
    # Convert back to string just to be safe
    genes_unique$chromosome_name <- as.character(genes_unique$chromosome_name)

  }  else {
    stop("No other option then 'first' currenly implementedq")
  }
  
  # Write to output file
  write.table(genes_unique, file = paste0(output_file, "_ensembl.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = c(attr, "final_gene_name"))
  
  # Make the gene linkers
  gene_linker <- genes_unique[,c("final_gene_name", "ensembl_gene_id")]
  colnames(gene_linker) <- c("hgnc_symbol", "ensembl_gene_id")

  write.table(gene_linker, file = paste0(output_file, "_name_to_ensembl.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE)
              
  write.table(gene_linker[,c(2,1)], file = paste0(output_file, "_ensembl_to_name.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE)
              
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
  #merging_strat = args[3]
  merging_strat="first"
  fetch_unique_genes(version, outfile, merging_strat)
}
