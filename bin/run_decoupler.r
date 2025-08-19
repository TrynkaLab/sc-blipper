#!/usr/bin/env Rscript

library(data.table)
library(decoupleR)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(optparse)
library(tibble)
library(tidyr)
library(OmnipathR)

#------------------------------------------------------------------------------
#' Simple heatmap with auto labels
#'
#' @param data A matrix
#' @param cellsize The size of the cells in x and y
#' @param cellwidth The size of the cells in x
#' @param cellheight THe size of the cells in y
#' @param limits Clip the data to upper and lower, must be a numeric of length 2
#' @param cluster Should clustering be applied
#' @param range Automatic setting for color range: 'symmetric', 'absolute', 'auto'
#' @param palette Palette to use in \code{\link{RColorBrewer::brewer.pal()}}
#' @param border Border color
#' @param ... Remaining arguments passed to  \code{\link{pheatmap::pheatmap()}}
#'
#' @returns The pheatmap plot
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @rdname heatmaps 
#' @export
plot_simple_hm <- function(data, cellsize = -1, cellwidth = 12, cellheight = 12, limits = NULL, cluster = T, range = "symmetric", palette = NULL, border = NA, silent=T, ...) {
  if (cellsize > 0) {
    cellwidth <- cellsize
    cellheight <- cellsize
  }
  
  if (is.null(limits)) {
    max <- max(data, na.rm=T)
    min <- min(data, na.rm=T)
    max.abs <- max(abs(data), na.rm=T)
    min.abs <- min(abs(data), na.rm=T)
  } else {
    if (length(limits) == 2 && is(limits, "numeric")) {
      max <- limits[2]
      min <- limits[1]
      max.abs <- max(abs(limits), na.rm=T)
      min.abs <- min(abs(limits), na.rm=T)
      range <- "auto"
    } else {
      stop("Limits must be a vector of length 2")
    }
  }
  
  if (range == "symmetric") {
    break.list <- seq(-max.abs, max.abs, by = max.abs / 50)
    if (is.null(palette)) {
      palette <- "RdBu"
    }
    cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
  } else if (range == "absolute") {
    if (is.null(palette)) {
      palette <- "Reds"
    }
    break.list <- seq(min, max.abs, by = (max.abs - min.abs) / 100)
    cols <- grDevices::colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
  } else if (range == "auto") {
    break.list <- seq(min, max, by = (max - min) / 100)
    if (is.null(palette)) {
      palette <- "RdBu"
    }
    cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
  } else {
    stop("Range must be 'symmetric', 'auto', or 'absolute'\n")
  }
  
  if (!cluster) {
    res <- pheatmap::pheatmap(data,
                              breaks = break.list,
                              col = cols,
                              cellwidth = cellwidth,
                              cellheight = cellheight,
                              border = border,
                              cluster_rows = F,
                              cluster_cols = F,
                              silent = silent,
                              ...
    )
  } else {
    res <- pheatmap::pheatmap(data,
                              breaks = break.list,
                              col = cols,
                              cellwidth = cellwidth,
                              cellheight = cellheight,
                              border = border,
                              silent = silent,
                              ...
    )
  }
  
  return(ggplotify::as.ggplot(res$gtable))
}



#-------------------------------------------------
# Define command-line options
option_list <- list(
  make_option(c("-m", "--matrix"),
              type = "character", help = "Path to numeric matrix CSV file (genes in rows by default)."),
  make_option(c("-o", "--output_prefix"),
              type = "character", help = "Prefix for output files."),
  make_option(c("-t", "--transpose"),
              type = "logical", default = FALSE,
              help = "Whether to transpose the matrix before processing. Default: FALSE."),
  make_option(c("--cache_dir"),
              type = "character", default = NULL,
              help = "Cache dir for OmniPathR. Default: NULL (home dir)"),
  make_option(c("--id_linker"),
            type = "character", default = NULL,
            help = "TSV file with gene name to ensembl (old, new) to convert ominpath with Default: NULL")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# Set the cache folder for OmniPathR
if (!is.null(opt$cache_dir)) {
  OmnipathR::omnipath_set_cachedir(opt$cache_dir)
}

# Check required arguments
if (any(is.null(opt$matrix), is.null(opt$output_prefix))) {
  print_help(opt_parser)
  stop("Error: --matrix, and --output_prefix are required.", call. = FALSE)
}

# If the input has ensembl ids use this file to convert omnipath
if (!is.null(opt$id_linker)) {
  mapping               <- read.table(opt$id_linker, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(mapping)     <- c("old", "new")
  gene_replacement      <- setNames(mapping$new, mapping$old)
  gene_replacement_rev  <- setNames(mapping$old, mapping$new)
} else {
  gene_replacement  <- NULL
  gene_replacement_rev <- NULL
}

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Get progeny information
net <- get_progeny(organism = 'human', top = 500)

if (!is.null(gene_replacement)) {
  net[,"target"] <- gene_replacement[net[,"target"]]
}

mat          <- na.omit(mat)

# Infer activities
activities <- run_mlm(mat=mat,
                      net=net,
                      .source='source',
                      .target='target',
                      .mor='weight',
                      minsize = 5)

activity.matrix <- activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')


write.table(activities, file=paste0(opt$output_prefix, "_progeny_activities.tsv"), quote=F, row.names=F, sep="\t")
write.table(activity.matrix, file=paste0(opt$output_prefix, "_progeny_activity_matrix.tsv"), quote=F, row.names=F)

#-------------------------------------------------------------------------------
# Visualization
activity.matrix.p <- activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'p_value') %>%
  column_to_rownames('source')

# Select the significant pathways (for progeny there are only a few, so just plot em all)
#sig.tfs <- rowSums(activity.matrix.p < 0.05/nrow(activities)) >=1 
#df.plot <- activity.matrix[sig.tfs, ]
df.plot <- activity.matrix

# Plotting ensembl ids is not very usefull, so in case they are , convert back to gene names
if (!is.null(gene_replacement_rev)) {
  rownames(df.plot) <- gene_replacement_rev[rownames(df.plot)]
}

pdf(width=2+(nrow(df.plot)*0.5), height=2+(ncol(df.plot)*0.5), file=paste0(opt$output_prefix, "_pathways.pdf"))
plot_simple_hm(df.plot, cellheight = 15, cellwidth = 15)
dev.off()


pdf(width=2+(nrow(df.plot)*0.5), height=2+(ncol(df.plot)*0.5), file=paste0(opt$output_prefix, "_pathways_scaled.pdf"))
plot_simple_hm(t(scale(t(df.plot))), cellheight = 15, cellwidth = 15)
dev.off()

#-------------------------------------------------------------------------------
# Get CollectTRI information
net <- get_collectri(organism='human', split_complexes=FALSE)

if (!is.null(gene_replacement)) {
  net[,"target"] <- gene_replacement[net[,"target"]]
}

mat          <- na.omit(mat)

# Infer activities
activities <- run_ulm(mat=mat,
                      net=net,
                      .source='source',
                      .target='target',
                      .mor='mor',
                      minsize = 5)

activity.matrix <- activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')

write.table(activities, file=paste0(opt$output_prefix, "_collectri_activities.tsv"), quote=F, row.names=F, sep="\t")
write.table(activity.matrix, file=paste0(opt$output_prefix, "_collectri_activity_matrix.tsv"), quote=F, row.names=F)

#-------------------------------------------------------------------------------
# Visualization
activity.matrix.p <- activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'p_value') %>%
  column_to_rownames('source')

# Select significant ones
sig.tfs <- rowSums(activity.matrix.p < 0.05/nrow(activities)) >=1 
df.plot <- activity.matrix[sig.tfs, ]

# Plotting ensembl ids is not very usefull, so in case they are , convert back to gene names
if (!is.null(gene_replacement_rev)) {
  rownames(df.plot) <- gene_replacement_rev[rownames(df.plot)]
}

pdf(width=2+(nrow(df.plot)*0.2), height=2+(ncol(df.plot)*1), file=paste0(opt$output_prefix, "_signif_tfs.pdf"))
plot_simple_hm(df.plot, cellheight = 8, cellwidth = 20)
dev.off()

pdf(width=2+(nrow(df.plot)*0.2), height=2+(ncol(df.plot)*1), file=paste0(opt$output_prefix, "_signif_tfs_scaled.pdf"))
plot_simple_hm(t(scale(t(df.plot))), cellheight = 8, cellwidth = 20)
dev.off()