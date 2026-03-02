#!/usr/bin/env Rscript
# summarize_cnmf_qc.r
#
# Summarize CNMF QC/annotation files and plot clustering metrics.
#
# Usage:
#   Rscript summarize_cnmf_qc.r "file1.tsv,file2.tsv,..." output_prefix
#
# Arguments:
#   1) A comma-separated list of input files (full paths). Example:
#        "/path/to/a.annotation.tsv,/path/to/b.annotation.tsv"
#      (If you prefer, you can build the comma-separated list using shell expansion.)
#   2) Output prefix: files will be written as <prefix>.pdf and <prefix>.annotation.combined.tsv
#
# The script reads all annotation tables, combines them by row, aggregates
# a set of QC metrics by cluster number `k` using `min()`, and writes a
# PDF plot and a combined annotation TSV.

library(data.table)
library(ggplot2)
library(optparse)

#---------------------------------------------------------
# Plot clustering metrics in a faceted plot
plot_clustering_metrics <- function(df) {
  # Select numeric columns excluding 'k'
  numeric_cols <- sapply(df, is.numeric)
  numeric_cols[names(df) == "k"] <- FALSE
  metric_cols <- names(df)[numeric_cols]
  
  # Convert to long format using base R reshape
  df_long <- reshape(df[, c("k", metric_cols)], 
                     direction = "long", 
                     varying = metric_cols,
                     v.names = "value",
                     times = metric_cols,
                     timevar = "metric")
  
  # Remove row names (they become time column index)
  rownames(df_long) <- NULL
  df_long$metric <- factor(df_long$metric, levels=metric_cols)
  
  # Create faceted plot
  p1 <- ggplot(df_long, aes(x = k, y = value)) +
    geom_point(size = 1, color = "darkgrey") +
    geom_line(color = "darkgrey", linewidth = 0.5) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    labs(
      x = "Number of Clusters (k)",
      y = "Metric Value",
      title = "Clustering Metrics by Number of Clusters"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),                 # no full box
      axis.line.x = element_line(color = "black"),    # x-axis bar
      axis.line.y = element_line(color = "black"),    # y-axis bar
      axis.ticks = element_line(color = "black")      # axis tick bars
    )
    
  return(p1)
}

# Helper: read table and optionally set first column as rownames
read <- function(file, rn = TRUE) {
  x <- fread(file, data.table = FALSE)
  if (rn && ncol(x) >= 1) {
    rownames(x) <- x[, 1]
    if (ncol(x) > 1) x <- x[, -1, drop = FALSE]
    else x <- data.frame(row.names = rownames(x))
  }
  return(x)
}
#---------------------------------------------------------
# Use optparse for command-line argument parsing
option_list <- list(
  make_option(c("-i", "--inputs"), type = "character",
              help = "Comma-separated list of input annotation files (full paths)",
              metavar = "FILES"),
  make_option(c("-o", "--out"), type = "character",
              help = "Output prefix for produced files (no extension)",
              metavar = "PREFIX")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$inputs) || is.null(opt$out)) {
  print_help(parser)
  stop("Both --inputs and --out are required.")
}

input_arg <- opt$inputs
out_prefix <- opt$out

#---------------------------------------------------------
# split comma-separated list and trim whitespace
files <- strsplit(input_arg, ",")[[1]]
files <- trimws(files)

# verify existence of files
exists_mask <- file.exists(files)
if (!all(exists_mask)) {
  warning("Some input files were not found:\n", paste(files[!exists_mask], collapse = ", "))
  files <- files[exists_mask]
}
if (length(files) == 0) {
  stop("No valid input files provided.")
}

# Read and combine annotation tables
annots <- list()
for (file in files) {
  annots[[basename(file)]] <- read(file)
}
annot <- do.call(rbind, annots)

# Aggregate metrics by cluster number `k` using the minimum as before
metric_cols      <- c("run_r2", "run_sse", "run_iter_count_perc", "run_silhouette", "run_calinski_harabasz", "run_davies_bouldin", "run_mean_density", "run_median_density", "min_edist")
present_metrics  <- intersect(metric_cols, colnames(annot))
if (length(present_metrics) == 0) stop("No expected metric columns found in annotations.")

df.plot <- aggregate(annot[, present_metrics, drop = FALSE], by = list(k = annot[, "k"]), FUN = min)

# Write combined annotation table
combined_out <- paste0(out_prefix, "_k_summary.tsv")
write.table(df.plot, file = combined_out, sep = "\t", quote = FALSE, col.names = NA)


# Create PDF of the metrics plot
pdf(file = paste0(out_prefix, "_k_summary.pdf"), width = 10, height = 5)
print(plot_clustering_metrics(df.plot))
dev.off()

message("Wrote: ", combined_out)
message("Wrote: ", paste0(out_prefix, "_k_summary.pdf"))

