#!/usr/bin/env Rscript

library(data.table)
library(decoupleR)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(optparse)

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
              help = "Whether to transpose the matrix before processing. Default: FALSE.")
)


# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# Check required arguments
if (any(is.null(opt$matrix), is.null(opt$output_prefix))) {
  print_help(opt_parser)
  stop("Error: --matrix, and --output_prefix are required.", call. = FALSE)
}

#-------------------------------------------------------------------------------
# Load numeric matrix (genes in rows, columns are samples or contrasts)
mat           <- as.matrix(fread(opt$matrix, data.table=FALSE))
rownames(mat) <- mat[,1]
mat           <- mat[,-1]

# Transpose if requested
if (opt$transpose) {
  mat <- t(mat)
}

#-------------------------------------------------------------------------------
# Get progeny information
net <- get_progeny(organism = 'human', top = 500)

# Select only overlapping genes and convert to matrix
unique.genes <- unlist(unique(net[,"target"]))
ol           <- intersect(unique.genes, rownames(mat))
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

# Select the significant pathways
#sig.tfs <- rowSums(activity.matrix.p < 0.05/nrow(activities)) >=1 
#df.plot <- activity.matrix[sig.tfs, ]
df.plot <- activity.matrix

# Scale per pathway
#df.plot <- t(scale(t(df.plot)))

pdf(width=10, height=5, file=paste0(opt$output_prefix, "_pathways.pdf"))
plot_simple_hm(df.plot, cellheight = 15, cellwidth = 15)
dev.off()


pdf(width=10, height=5, file=paste0(opt$output_prefix, "_pathways_scaled.pdf"))
plot_simple_hm(t(scale(t(df.plot))), cellheight = 15, cellwidth = 15)
dev.off()


#-------------------------------------------------------------------------------
# Get CollectTRI information
net <- get_collectri(organism='human', split_complexes=FALSE)

# Select only overlapping genes and convert to matrix
unique.genes <- unlist(unique(net[,"target"]))
mat          <- na.omit(mat)
ol           <- intersect(unique.genes, rownames(mat))

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

write.table(activities, file=paste0(opt$output_prefix, "_collectri_activities.tsv"), quote=F, row.names=T, sep="\t")
write.table(activity.matrix, file=paste0(opt$output_prefix, "_collectri_activity_matrix.tsv"), quote=F, row.names=T)

#-------------------------------------------------------------------------------
# Visualization
activity.matrix.p <- activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'p_value') %>%
  column_to_rownames('source')

# Select significant ones
sig.tfs <- rowSums(activity.matrix.p < 0.05/nrow(activities)) >=1 
df.plot <- activity.matrix[sig.tfs, ]

pdf(width=10, height=20, file=paste0(opt$output_prefix, "_signif_tfs.pdf"))
plot_simple_hm(df.plot, cellheight = 2, cellwidth = 20)
dev.off()

pdf(width=10, height=20, file=paste0(opt$output_prefix, "_signif_tfs_scaled.pdf"))
plot_simple_hm(t(scale(t(df.plot))), cellheight = 2, cellwidth = 20)
dev.off()