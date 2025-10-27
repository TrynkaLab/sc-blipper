#! /usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)


# Read gnomAD constraint metrics
gnomad_lof <- read_tsv("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv")
gnomad_lof <- gnomad_lof[gnomad_lof$canonical,]
gnomad_lof <- gnomad_lof[startsWith(gnomad_lof$gene_id, "ENSG"),]

# Function to write GMT file
write_gmt <- function(gene_list, description, output_file) {
  write(paste(description, description, paste(gene_list, collapse="\t")), 
        file=output_file, 
        append=TRUE)
}

# Process data
genes_df <- gnomad_lof %>%
  filter(!is.na(lof.pLI)) %>%
  arrange(desc(lof.pLI)) %>%
  select(gene, gene_id, lof.pLI)

# Initialize output files
file_path <- "gene_sets/"

if (file.exists(paste0(file_path, "symbols/gnomad.v4.lof.symbols.gmt"))) {
  file.remove(paste0(file_path, "symbols/gnomad.v4.lof.symbols.gmt"))
}

if (file.exists(paste0(file_path, "ensembl/gnomad.v4.lof.ensembl.gmt"))) {
  file.remove(paste0(file_path, "ensembl/gnomad.v4.lof.ensembl.gmt"))
}

# Write different lof.pLI threshold gene sets
lof.pLI_thresholds <- c(0.8, 0.9, 0.95, 0.99)
for (thresh in lof.pLI_thresholds) {
  genes_thresh <- genes_df %>%
    filter(lof.pLI >= thresh)
  
  # Write version with ENSEMBL IDs
  write_gmt(genes_thresh$gene_id,
            sprintf("lof.pLI_%.2f", thresh),
            paste0(file_path, "ensembl/gnomad.v4.lof.ensembl.gmt"))
  
  # Write version with gene symbols
  write_gmt(genes_thresh$gene,
            sprintf("lof.pLI_%.2f", thresh),
            paste0(file_path, "symbols/gnomad.v4.lof.symbols.gmt"))
}

# Write top N genes by lof.Z
# Process data
genes_df <- gnomad_lof %>%
  filter(!is.na(lof.z_score)) %>%
  arrange(desc(lof.z_score)) %>%
  select(gene, gene_id, lof.z_score)

top_n_values <- c(50, 100, 250, 500)
for (n in top_n_values) {
  top_genes <- genes_df %>%
    slice_head(n = n)
  
  # Write version with ENSEMBL IDs
  write_gmt(top_genes$gene_id,
            sprintf("lof.z_score_top_%d", n),
            paste0(file_path, "ensembl/gnomad.v4.lof.ensembl.gmt"))
  
  # Write version with gene symbols
  write_gmt(top_genes$gene,
            sprintf("lof.z_score_top_%d", n),
            paste0(file_path, "symbols/gnomad.v4.lof.symbols.gmt"))
}

