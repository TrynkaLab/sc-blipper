#!/usr/bin/env Rscript

library(Seurat)
library(anndata)
library(reticulate)

# This is where the reticulate python lives
Sys.setenv(RETICULATE_PYTHON=Sys.which("python"))

convert_cols_to_string_or_number <- function(df) {
  df[] <- lapply(df, function(col) {
    num_col <- suppressWarnings(as.numeric(as.character(col)))
    if (sum(is.na(num_col)) > 0.5 * length(num_col)) {
      as.character(col)
    } else {
      num_col
    }
  })
  return(df)
}

# Function to update gene IDs using a mapping file and save unupdated IDs
update_gene_ids <- function(counts_mat, mapping_file, unupdated_file = "unmapped_genes.txt") {
  mapping           <- read.table(mapping_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(mapping) <- c("old", "new")
  gene_replacement  <- setNames(mapping$new, mapping$old)
  genes             <- rownames(counts_mat)
  
  # Identify genes that have replacements
  has_replacement <- genes %in% names(gene_replacement)
  
  # Collect unupdated gene IDs
  unupdated_genes <- genes[!has_replacement]
  
  # Write unupdated gene IDs to file
  if (length(unupdated_genes) > 0) {
    write.table(unupdated_genes, file=unupdated_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
    warning(length(unupdated_genes), " gene IDs did not get updated and were saved to ", unupdated_file)
  } else {
    message("All gene IDs were updated.")
  }

  # Perform the replacement
  genes_updated        <- ifelse(has_replacement, gene_replacement[genes], genes)
  rownames(counts_mat) <- genes_updated
  return(counts_mat)
}


convert_seurat_to_h5ad_counts <- function(input_rds, output_h5ad, mapping_file=NULL, unupdated_file=NULL) {
    seu        <- LoadSeuratRds(input_rds)
    counts_mat <- Seurat::GetAssayData(seu, assay = DefaultAssay(seu), layer = "counts")

    # Save the original gene ids
    gene.ids.orig <- rownames(counts_mat)
    
    # Apply gene ID mapping if provided
    if (!is.null(mapping_file)) {
      counts_mat <- update_gene_ids(counts_mat, mapping_file, unupdated_file)
    }

    cell_meta <- seu@meta.data
    cell_meta <- convert_cols_to_string_or_number(cell_meta)
    if (!inherits(counts_mat, "dgCMatrix")) {
        counts_mat <- as(counts_mat, "dgCMatrix")
    }

    ad <- AnnData(
        X = counts_mat,
        var = data.frame(cell_meta),
        obs = data.frame(row.names=rownames(counts_mat), id=rownames(counts_mat), orig=gene.ids.orig)
    )
    
    # Transpose to comply with scanpy standard, this is MUCH faster in python
    ad <- ad$transpose()
    ad$write_h5ad(output_h5ad)
    message("Done! Saved counts + cell metadata to: ", output_h5ad)
}

# Command line interface
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2 || length(args) > 3) {
        cat("Usage: Rscript seurat_to_h5ad_counts.R <input_seurat_rds> <output_h5ad> [gene_id_mapping.txt [unupdated_file.txt]]\n")
        quit(status = 1)
    }
    mapping_file <- if (length(args) == 3) args[3] else NULL
    
    unupdated_file <- if (length(args) == 4) args[4] else NULL
    convert_seurat_to_h5ad_counts(args[1], args[2], mapping_file, unupdated_file)
}

