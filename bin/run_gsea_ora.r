#!/usr/bin/env Rscript

library(fgsea)
library(data.table)
library(tools)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-m", "--matrix"),
              type = "character", help = "Path to numeric matrix CSV file (genes in rows by default)"),
  make_option(c("-g", "--gmt"),
              type = "character", help = "Comma-separated list of GMT files."),
  make_option(c("-o", "--output_prefix"),
              type = "character", help = "Prefix for output files."),
  make_option(c("-t", "--transpose"),
              type = "logical", default = FALSE,
              help = "Whether to transpose the matrix before processing. Default: FALSE."),
  make_option(c("--threshold"),
              type = "numeric", default = 0.05,
              help = "Threshold to binarize the input matrix. Use pvalue/fdr threshold if -m is pvalues, use 1 if input is binary, use 0 if input is cNMF facotors. Interacts with --absolute Default: 0.05"),
  make_option(c("--use_top"),
              type = "character", default = NULL,
              help = "Instead of binzarizing, rank each collumn and use the top x genes (absolute if --absolute is specified). Overrides --threshold. Comma seperated list '100,200,500'. Default: NULL"),
  make_option(c("--absolute"),
              type = "boolean", default = TRUE,
              help = "When --threshold or --use_top, use the absolute value instead. Default: TRUE"),
  make_option(c("-u", "--universe"),
              type = "character", default = NA,
              help = "Optional gene universe file (one gene per line). Defaults to all genes in the matrix.")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (!is.null(opt$use_top)) {
  opt$use_top <- as.numeric(strsplit(opt$use_top, split=",")[[1]])
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

results  <- list()

for (gmt_file in gmt_files) {
  pathways <- prepare_pathways(gmt_file, gene_universe)
  gmt_name <- gsub(".symbols.gmt$|.gmt$|.ensembl.gmt$", "", basename(gmt_file))
  
  for (col_name in colnames(mat)) {
    cur.mat <- mat[, col_name]
    
    if (opt$absolute) {
      sign    <- sign(cur.mat)
      cur.mat <- abs(cur.mat)
    }
    
    if (!is.null(opt$use_top)) {
      rank <- sort(cur.mat, decreasing=T)
      
      for (ngene in opt$use_top) {
        
        if (ngene >= length(rank)) {
            warning("More genes requested then are in the genelist. Skipping this iteration")
            next()
        }
        
        fgsea_res <- fora(pathways = pathways,
                          genes = names(rank)[1:ngene],
                          universe=gene_universe,
                          minSize = 5,
                          maxSize = 500) 
        fgsea_res$condition  <- col_name
        fgsea_res$condition2 <- paste0("top",ngene)
        fgsea_res$database   <- gmt_name
        results[[paste0(col_name, "_", gmt_name, "_",paste0("top",ngene))]]  <- fgsea_res
      }
      
    } else {
      genes <- sign[cur.mat > opt$threshold] * cur.mat[cur.mat > opt$threshold]
      
      # ALl genes, up and down
      if (length(genes) > 0 ) {
        fgsea_res <- fora(pathways = pathways,
                          genes = names(genes),
                          universe=gene_universe,
                          minSize = 5,
                          maxSize = 500) 
        fgsea_res$condition  <- col_name
        fgsea_res$condition2 <- "ALL"
        fgsea_res$database   <- gmt_name
        results[[paste0(col_name, "_", gmt_name, "_ALL")]]  <- fgsea_res
      }
      
      # UP genes
      cur.genes <- names(genes[genes >0])
      if (length(cur.genes) > 0 ) {
        fgsea_res <- fora(pathways = pathways,
                          genes = cur.genes,
                          universe=gene_universe,
                          minSize = 5,
                          maxSize = 500) 
        fgsea_res$condition  <- col_name
        fgsea_res$condition2 <- "UP"
        fgsea_res$database   <- gmt_name
        results[[paste0(col_name, "_", gmt_name, "_UP")]]  <- fgsea_res
      }
      
      # DOWN genes
      cur.genes <- names(genes[genes < 0])
      
      if (length(cur.genes) > 0 ) {
        fgsea_res <- fora(pathways = pathways,
                          genes = cur.genes,
                          universe=gene_universe,
                          minSize = 5,
                          maxSize = 500) 
        fgsea_res$condition  <- col_name
        fgsea_res$condition2 <- "DOWN"
        fgsea_res$database   <- gmt_name
        results[[paste0(col_name, "_", gmt_name, "_DOWN")]]  <- fgsea_res
      }
    }
  }
  
}

res                  <- as.data.frame(do.call(rbind, results))
res[,"overlapGenes"] <- sapply(res[,"overlapGenes"], paste0, collapse=",")
out_file             <- paste0(opt$output_prefix, "_ora_results.tsv")
write.table(res, out_file, sep="\t", quote=F, row.names=F)

