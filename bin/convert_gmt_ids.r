#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript convert_gmt_ids.r <input_gmt> <mapping_tsv> <output_gmt>\n")
  quit(save = "no", status = 1)
}

gmt_file     <- args[1]
mapping_file <- args[2]
output_file  <- args[3]

# -----------------------------------------------------------
read_gmt <- function(file) {
  lines <- readLines(file)
  strsplit(lines, "\t", fixed = TRUE)
}

write_gmt <- function(gmt_list, file) {
  lines <- vapply(gmt_list, function(x) paste(x, collapse = "\t"), character(1))
  writeLines(lines, file)
}

# -----------------------------------------------------------
# Load GMT data
gmt_data <- read_gmt(gmt_file)

# Load mapping table (header assumed: GeneSymbol    EnsemblID)
mapping <- read.table(mapping_file, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, quote = "")

# Create lookup: names = symbols, values = Ensembl IDs
symbol2ensembl <- setNames(mapping[,2], mapping[,1])

# Convert each pathway's genes
gmt_ensembl <- lapply(gmt_data, function(path) {
  pathway_info <- path[1:2]
  genes_symbols <- path[-c(1, 2)]
  
  # Map to Ensembl IDs and drop NAs
  genes_ensembl <- symbol2ensembl[genes_symbols]
  genes_ensembl <- unique(genes_ensembl[!is.na(genes_ensembl)])
  
  c(pathway_info, genes_ensembl)
})

# Write the converted GMT
write_gmt(gmt_ensembl, output_file)

cat("New GMT file with updated IDs saved to:", output_file, "\n")
