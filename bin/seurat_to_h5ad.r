#!/usr/bin/env Rscript

library(Seurat)
library(anndata)
library(reticulate)


# This is where the reticulate python lives, this should work with 
# conda, softpack modules and singularity containers
Sys.setenv(RETICULATE_PYTHON=Sys.which("python"))

# This is needed as h5ad only accepts prmitice types
convert_cols_to_string_or_number <- function(df) {
  df[] <- lapply(df, function(col) {
    # Try converting to numeric
    num_col <- suppressWarnings(as.numeric(as.character(col)))
    # If conversion results in too many NAs (non-numeric), keep as character
    if (sum(is.na(num_col)) > 0.5 * length(num_col)) {
      # Return as character
      as.character(col)
    } else {
      # Return as numeric
      num_col
    }
  })
  return(df)
}

# Convert and save the h5 file
convert_seurat_to_h5ad_counts <- function(input_rds, output_h5ad) {
    # Load Seurat object
    seu <- LoadSeuratRds(input_rds)
    
    # Extract counts matrix from default assay as sparse matrix
    counts_mat <- Seurat::GetAssayData(seu, assay = DefaultAssay(seu), layer = "counts")
    
    # Extract cell metadata (colData)
    cell_meta <- seu@meta.data
    
    # Convert to numeric or string primitives
    cell_meta <- convert_cols_to_string_or_number(cell_meta)
    
    # Convert counts to dgCMatrix (which anndata can handle through reticulate/Scipy)
    if (!inherits(counts_mat, "dgCMatrix")) {
        counts_mat <- as(counts_mat, "dgCMatrix")
    }
    
    # Create AnnData object
    ad <- AnnData(
        X = counts_mat,
        var = data.frame(cell_meta),
        obs = data.frame(row.names=rownames(counts_mat), id=rownames(counts_mat)),
    )
    
    # Transpose to follow anndata convention
    ad <- ad$transpose()

    # Save to h5ad
    ad$write_h5ad(output_h5ad)
    message("✅ Done! Saved counts + cell metadata to: ", output_h5ad)
        
}

# Example usage:
# convert_seurat_to_h5ad_counts("mydata.rds", "mydata_counts_only.h5ad")

# Command line interface
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) {
        cat("Usage: Rscript seurat_to_h5ad_counts.R <input_seurat_rds> <output_h5ad>\n")
        quit(status = 1)
    }
    convert_seurat_to_h5ad_counts(args[1], args[2])
}
