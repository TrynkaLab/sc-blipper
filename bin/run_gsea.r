#!/usr/bin/env Rscript

library(fgsea)
library(data.table)
library(tools)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-m", "--matrix"),
              type = "character", help = "Path to numeric matrix CSV file (genes in rows by default)."),
  make_option(c("-g", "--gmt"),
              type = "character", help = "Comma-separated list of GMT files."),
  make_option(c("-o", "--output_prefix"),
              type = "character", help = "Prefix for output files."),
  make_option(c("-t", "--transpose"),
              type = "logical", default = FALSE,
              help = "Whether to transpose the matrix before processing. Default: FALSE."),
  make_option(c("-u", "--universe"),
              type = "character", default = NA,
              help = "Optional gene universe file (one gene per line). Defaults to all genes in the matrix.")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# Check required arguments
if (any(is.null(opt$matrix), is.null(opt$gmt), is.null(opt$output_prefix))) {
  print_help(opt_parser)
  stop("Error: --matrix, --gmt, and --output_prefix are required.", call. = FALSE)
}

# Load numeric matrix (genes in rows, columns are samples or contrasts)
mat           <- as.matrix(fread(opt$matrix, data.table=FALSE))
rownames(mat) <- mat[,1]
mat           <- mat[,-1]

# Transpose if requested
if (opt$transpose) {
  mat <- t(mat)
}

# Determine gene universe
if (!is.na(opt$universe)) {
  gene_universe <- scan(opt$universe, what = character())
  gene_universe <- intersect(rownames(mat), gene_universe)
} else {
  gene_universe <- rownames(mat)
}

# Function to read and filter pathways
prepare_pathways <- function(gmt_file, gene_universe) {
  pathways          <- gmtPathways(gmt_file)
  pathways_filtered <- lapply(pathways, function(genes) intersect(genes, gene_universe))
  pathways_filtered <- Filter(function(x) length(x) >= 5, pathways_filtered)
  return(pathways_filtered)
}

# Process each GMT file
gmt_files <- strsplit(opt$gmt, ",")[[1]]

for (gmt_file in gmt_files) {
  pathways <- prepare_pathways(gmt_file, gene_universe)
  gmt_name <- gsub(".symbols.gmt$|.gmt$|.ensembl.gmt$", "", basename(gmt_file))
  results  <- list()
  
  for (col_name in colnames(mat)) {
    ranks <- mat[, col_name]
    names(ranks) <- rownames(mat)
    ranks <- ranks[intersect(names(ranks), gene_universe)]
    ranks <- sort(ranks, decreasing = TRUE)
    
    fgsea_res <- fgsea(pathways = pathways,
                       stats = ranks,
                       minSize = 5,
                       maxSize = 500)
    
    fgsea_res$condition  <- col_name
    fgsea_res$condition2 <- "GSEA"
    results[[col_name]]  <- fgsea_res
  }
  
  res                 <- as.data.frame(do.call(rbind, results))
  res[,"leadingEdge"] <- sapply(res[,"leadingEdge"], paste0, collapse=",")
  out_file            <- paste0(opt$output_prefix, "_", gmt_name, "_fgsea_results.tsv")
  write.table(res, out_file, sep="\t", quote=F, row.names=F)
}

