#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(anndata)
  library(reticulate)
})

# Ensure correct Python environment
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
  df
}

# Update gene IDs using mapping file
update_gene_ids <- function(counts_mat, mapping_file, unupdated_file = "unmapped_genes.txt") {
  mapping           <- read.table(mapping_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(mapping) <- c("old", "new")
  gene_replacement  <- setNames(mapping$new, mapping$old)
  genes             <- rownames(counts_mat)
  
  has_replacement  <- genes %in% names(gene_replacement)
  unupdated_genes  <- genes[!has_replacement]
  
  if (length(unupdated_genes) > 0) {
    write.table(unupdated_genes, file=unupdated_file, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    warning(length(unupdated_genes), " gene IDs did not get updated and were saved to ", unupdated_file)
  } else {
    message("All gene IDs were updated.")
  }

  rownames(counts_mat) <- ifelse(has_replacement, gene_replacement[genes], genes)
  return(list(counts=counts_mat, mapping=gene_replacement))
}

# Filter genes based on supplied list (after ID mapping)
filter_genes <- function(counts_mat, filter_file, gene_mapping=NULL) {
  filter_genes <- read.table(filter_file, header=FALSE, stringsAsFactors=FALSE)[,1]
  
  # If mapping provided, update filter gene IDs too
  if (!is.null(gene_mapping)) {
    filter_genes <- ifelse(filter_genes %in% names(gene_mapping),
                           gene_mapping[filter_genes],
                           filter_genes)
  }
  
  keep <- rownames(counts_mat) %in% filter_genes
  return(counts_mat[keep, , drop=FALSE])
}

convert_seurat_to_h5ad_counts <- function(input_rds, output_h5ad, mapping_file=NULL, 
                                          unupdated_file=NULL, filter_file=NULL) {
  seu        <- LoadSeuratRds(input_rds)
  cell_meta  <- seu@meta.data
  counts_mat <- Seurat::GetAssayData(seu, assay = DefaultAssay(seu), layer = "counts")
  message("Read counts matrix of shape: ", dim(counts_mat))
  
  # Cleanup to save some memory
  rm(seu)
  gc()

  gene.ids.orig <- rownames(counts_mat)

  gene_mapping <- NULL
  if (!is.null(mapping_file)) {
    updated      <- update_gene_ids(counts_mat, mapping_file, unupdated_file)
    counts_mat   <- updated$counts
    gene_mapping <- updated$mapping
  }

  if (!is.null(filter_file)) {
    counts_mat <- filter_genes(counts_mat, filter_file, gene_mapping)
    message("Filtered counts down to: ", dim(counts_mat))
  }

  cell_meta <- convert_cols_to_string_or_number(cell_meta)
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }

  ad <- AnnData(
    X   = counts_mat,
    var = data.frame(cell_meta),
    obs = data.frame(row.names=rownames(counts_mat), id=rownames(counts_mat), orig=gene.ids.orig)
  )
  
  ad <- ad$transpose()
  ad$write_h5ad(output_h5ad)
  message("Done! Saved counts + cell metadata to: ", output_h5ad)
}

### Command line interface with optparse ###
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type="character", help="Input Seurat RDS file", metavar="file"),
    make_option(c("-o", "--output"), type="character", help="Output h5ad filename", metavar="file"),
    make_option(c("-m", "--mapping"), type="character", default=NULL, 
                help="Optional gene ID mapping file (tab-delimited)", metavar="file"),
    make_option(c("-u", "--unmapped"), type="character", default="unmapped_genes.txt", 
                help="File to save unmapped gene IDs [default: %default]", metavar="file"),
    make_option(c("-f", "--filter-genes"), type="character", default=NULL,
                help="Optional file with gene IDs to filter after mapping. One ID per line.", metavar="file")
  )

  parser <- OptionParser(option_list=option_list)
  opt <- parse_args(parser)

  if (is.null(opt$input) || is.null(opt$output)) {
    print_help(parser)
    quit(status=1)
  }

  convert_seurat_to_h5ad_counts(opt$input, opt$output, 
                                opt$mapping, opt$unmapped,
                                opt$`filter-genes`)
}
