#!/usr/bin/env Rscript

library(fgsea)
library(data.table)
library(tools)
library(optparse)
library(ggplot2)
#--------------------------------------------------
# Functions

# Take a matrix and return the top x values based on a grouping col
top_x_per_group <- function(mat, group_col, value_col, top_n = 5, absolute=FALSE, invert=FALSE) {
  df <- as.data.frame(mat)

  # If columns are numeric indices, convert to names
  if (is.numeric(group_col)) group_col <- names(df)[group_col]
  if (is.numeric(value_col)) value_col <- names(df)[value_col]

  # Ensure the value column is numeric
  df[[value_col]] <- as.numeric(df[[value_col]])
  
  # Split data.frame into a list by group
  split_list <- split(df, df[[group_col]])
  
  # For each group, sort by value and take top_n rows
  top_list <- lapply(split_list, function(subdf) {
    
    if (absolute) {
      subdf[order(-abs(subdf[[value_col]]), decreasing=invert), ][seq_len(min(top_n, nrow(subdf))), ]
    } else {
      subdf[order(-subdf[[value_col]], decreasing=invert), ][seq_len(min(top_n, nrow(subdf))), ]
    }
  })
  
  # Combine back into a single data.frame
  result <- do.call(rbind, top_list)
  
  rownames(result) <- NULL
  return(result)
}

# Function to read and filter pathways (in gmt format) to a universe
prepare_pathways <- function(gmt_file, gene_universe) {
  pathways          <- gmtPathways(gmt_file)
  pathways_filtered <- lapply(pathways, function(genes) intersect(genes, gene_universe))
  pathways_filtered <- Filter(function(x) length(x) >= 5, pathways_filtered)
  return(pathways_filtered)
}

cut_from_front <- function(strings, maxlen) {
  n <- nchar(strings)
  start_pos <- pmax(n - maxlen + 1, 1) # starting index
  substr(strings, start_pos, n)
}

dotplot <- function(df) {
  #df$Signif <- df$padj < 0.05  # create logical for padj significance
  df[,"Signif"] <- df[,"padj"] < 0.05
  
  p <- ggplot(
    df,
    aes(
      x = paste0(df$condition, "_", df$condition2),                       # or reorder(pathway, -foldEnrichment) for sorting
      y = cut_from_front(pathway, 40),                     # or another dimension, adjust as needed
      color = -log10(pval),
      shape = Signif
    )
  ) +
    geom_point(alpha = 0.7, size=3) +
    scale_color_gradient(low = "blue", high = "red", guide = guide_colorbar(title = "-log10 p-value")) +
    scale_shape_manual(values = c(16, 17), labels = c("padj >= 0.05", "padj < 0.05"), guide = guide_legend(title = "Significant")) +
    theme_bw() +
    labs(
      x = "Condition",
      y = "Pathway",
      title = "Dotplot of Pathways:, Color = p-value, Shape = padj < 0.05"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + facet_grid(database~., scales="free", space = "free")
    
    return(p)
    
    #    scale_size_continuous(range = c(2, 10), guide = guide_legend(title = "Fold Enrichment")) +

}

#-------------------------------------------------
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
              help = "Whether to transpose the matrix before processing so genes are on the rows. Default: FALSE."),
  make_option(c("-u", "--universe"),
              type = "character", default = NA,
              help = "Optional gene universe file (one gene per line). Defaults to all genes in the matrix."),
  make_option(c("--update_rows"),
              type = "character", help = "Path to file with 'old' 'new' ids. Aplied after transpose", default=NULL)
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
mat           <- fread(opt$matrix, data.table=FALSE, header=T)
rownames(mat) <- mat[,1]
mat           <- mat[,-1]
mat           <- as.matrix(mat)

if (!class(mat[1,1]) %in% c("numeric", "integer")) {
  stop("[ERROR] Input matrix is not numeric\n")
  q(save=FALSE, status=1)
}

# Transpose if requested
if (opt$transpose) {
  mat <- t(mat)
}

if (!is.null(opt$update_rows)) {
  mapping           <- read.table(mapping_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(mapping) <- c("old", "new")
  gene_replacement  <- setNames(mapping$new, mapping$old)
  rownames(mat)     <- gene_replacement[rownames(mat)]
}

# Determine gene universe
if (!is.na(opt$universe)) {
  gene_universe <- scan(opt$universe, what = character())
  gene_universe <- intersect(rownames(mat), gene_universe)
  cat("Universe-matrix overlap of ", length(gene_universe), "\n")
} else {
  gene_universe <- rownames(mat)
}

# Process each GMT file
gmt_files <- strsplit(opt$gmt, ",")[[1]]

results  <- list()

for (gmt_file in gmt_files) {
  pathways <- prepare_pathways(gmt_file, gene_universe)
  gmt_name <- gsub(".symbols.gmt$|.gmt$|.ensembl.gmt$", "", basename(gmt_file))
  
  for (col_name in colnames(mat)) {
    cat("debug - ", gmt_file, " - ", col_name, "\n")
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
    fgsea_res$database   <- gmt_name
    results[[paste0(col_name, "_", gmt_name)]]  <- fgsea_res
  }
}

res                 <- as.data.frame(do.call(rbind, results))
res[,"leadingEdge"] <- sapply(res[,"leadingEdge"], paste0, collapse=",")
out_file            <- paste0(opt$output_prefix, "_fgsea_results.tsv")
write.table(res, out_file, sep="\t", quote=F, row.names=F)


# Sort on top pvalues (always absolute, smallest first  )
res.sig       <- res
res.sig       <- res[res$padj < 0.05,]
res.sig$group <- paste0(res.sig$database, ";", res.sig$condition, ";", res.sig$condition2)
res.sig       <- res.sig[!is.na(res.sig$database),]
res.top       <- top_x_per_group(mat=res.sig,group_col="group", value_col="pval", absolute=T, invert=T)

out_file            <- paste0(opt$output_prefix, "_fgsea_results_top5.tsv")
write.table(res.top, out_file, sep="\t", quote=F, row.names=F)

pdf(width=(5+(0.5*length(unique(res.top$condition)))), height=(3 + (nrow(res.top)*0.1)), file=paste0(opt$output_prefix, "_fgsea_results_top5.pdf"))
dotplot(res.top)
dev.off()
