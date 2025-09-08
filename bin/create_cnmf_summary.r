#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(reshape2)
})

#-------------------------------------------------------------------------------
# Functions
read <- function(file, rn=T) {
  x <- fread(file, data.table=F)
  if (rn) {
    rownames(x) <- x[,1]
    x <- x[,-1, drop=F]
  }
  return(x)
}

get.topn <- function(df, sort.col="pvalue", group.col="group", top.n=5) {
  df[,sort.col] <- as.numeric(df[,sort.col])
  split_df <- split(df, df[[group.col]])
  top_n_list <- lapply(split_df, function(subdf) {
    subdf[order(-subdf[[sort.col]]), ][1:min(top.n, nrow(subdf)), ]
  })
  result_df <- do.call(rbind, top_n_list)
  return(result_df)
}

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


#-------------------------------------------------------------------------------
# Option Parsing
option_list <- list(
  make_option(c("-n", "--topn"), type="integer", default=5,
              help="Top N pathways to return [default: %default]"),

  make_option(c("-s", "--sep"), type="character", default=", ",
              help="Separator for gene names [default: '%default']"),

  make_option(c("-t", "--tests"), type="character", default=NULL,
              help="Comma-separated list of tests (e.g. 'MAGMA,ORA;top50,TF_TARGET')"),

  make_option(c("-d", "--databases"), type="character", default=NULL,
              help="Comma-separated list of databases"),

  make_option(c("-S", "--spectra"), type="character", default=NULL,
              help="Required: spectra file"),

  make_option(c("-e", "--enrichment"), type="character", default=NULL,
              help="Enrichment file or NULL"),
  
  make_option(c("--enrichmentThreshold"), type="numeric", default=5e-5,
            help="Pvalue threshold to consider enrichments [default: %default]"),

  make_option(c("-a", "--annot"), type="character", default=NULL,
              help="Annotation file or NULL"),

  make_option(c("--logSpectra"), action="store_true", default=FALSE,
              help="Log2-transform spectra matrix [default: %default]"),

  make_option(c("--scaleSpectra"), action="store_true", default=FALSE,
              help="Scale the spectra matrix [default: %default]"),

  make_option(c("-o", "--output"), type="character", default="output.tsv",
              help="Output prefix [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

#-------------------------------------------------------------------------------
# Process options into vectors/lists
tests <- if (!is.null(opt$tests)) unlist(strsplit(opt$tests, ",")) else NULL
databases <- if (!is.null(opt$databases)) unlist(strsplit(opt$databases, ",")) else NULL

spectra.file <- opt$spectra
enrichment.file <- opt$enrichment
annot.file <- opt$annot
top.n <- opt$topn
sep <- opt$sep
log.spectra <- opt$logSpectra
scale.spectra <- opt$scaleSpectra
output.prefix <- opt$output
enrich.threshold <- opt$enrichmentThreshold

if (is.null(spectra.file)) {
  stop("You must provide a spectra file with -S / --spectra")
}

#-------------------------------------------------------------------------------
# Spectra
spectra     <- read(spectra.file)
spectra     <- t(spectra)

# Find the top.n genes based on the spectra file
top.spectra <- apply(apply(spectra, 2, function(x){rownames(spectra)[order(x, decreasing = T)]})[1:top.n,], 2, paste0, collapse=sep)
top.mat     <- data.frame(top_gep_genes=top.spectra)
cat("[INFO] read spectra\n")
#-------------------------------------------------------------------------------
# Enrichment
if (!is.null(enrichment.file)) {
  enrich       <- read(enrichment.file, rn=F)
  cat("[INFO] read erichment\n")
  print(head(enrich))

  if (is.null(tests)) {
    tests <- unique(enrich$test)
  }
  
  if (is.null(databases)) {
    databases <- unique(enrich$database)
  }
  
  enrich       <- enrich[enrich$test %in% tests,]
  enrich       <- enrich[enrich$database %in% databases,]

  if (nrow(enrich) >=1) {
      enrich$group  <- paste0(enrich$test, "_", enrich$condition)
      enrich        <- enrich[,c("group", "trait", "condition", "test", "pvalue")]
      enrich        <- enrich[enrich$pvalue < enrich.threshold,]
      
      enrich$pvalue <- -log10(enrich$pvalue)
      
      top          <- get.topn(enrich, top.n=top.n)
      # Aggregate by group
      top.agg <- aggregate(top[,c("condition", "test", "trait")], by=list(top$group), function(x){paste0(unique(x), collapse = sep)})
      
      # Convert to matrix
      top.agg <- acast(top.agg, condition ~ test, value.var = "trait")
      
      top.mat <- cbind(top.mat, top.agg[colnames(spectra),])
  } else {
    warning("[WARN] No enrichment records left after filtering\n")
  } 

}

#-------------------------------------------------------------------------------
# Annotations
if (!is.null(annot.file)) {
  annot           <- read(annot.file)
    
  ol              <- intersect(rownames(spectra), rownames(annot))
  m               <- spectra[ol,]
  m2              <- aggregate(m, by=list(group=annot[ol,1]), FUN=mean)
  rownames(m2)    <- m2[,1]
  m2              <- t(m2[,-1])
  if (log.spectra) {
    m2 <- log2(m2+1)
  }
  if (scale.spectra) {
    m2 <- scale(m2)
  }
  
  # Cluster 
  d <- dist(t(m2))       # transpose to cluster columns
  hc <- hclust(d)   
  
  
  top.mat <- cbind(top.mat, m2[rownames(top.mat), hc$order])
  
  pdf(width=2+(nrow(m2)*0.5), height=2+(ncol(m2)*0.5), file=paste0(output.prefix, "_marker_heatmap.pdf"))
  plot(plot_simple_hm(m2))
  dev.off()
}

#-------------------------------------------------------------------------------
# Output
top.mat <- data.frame(rownames=rownames(top.mat), gep=paste0("GEP",rownames(top.mat)), top.mat)
write.table(top.mat,
            file=paste0(output.prefix, "_summary.tsv"),
            quote=F,
            sep="\t",
            row.names=F)
