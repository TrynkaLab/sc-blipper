#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(reshape2)
  library(ggplot2)
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

# Creates a profile plot for a given GEP
plot.top.k.genesets <- function(spectra, gep, k=20, cyto=NULL, tfs=NULL, scale="mad") {
  
  cur <- spectra[,gep]
  
  plots <- list()
  
  gepname <- paste0("GEP", gep)
  
  # Top k genes
  cur <- sort(cur, decreasing = T)
  genes <- names(cur)
  
  if (scale == "mad") {
    cur <- (cur - median(cur)) / mad(cur)
  } else if (scale == "z") {
    cur <- (cur - mean(cur)) / sd(cur)
  } else if (scale == "minmax") {
    cur <- (cur - min(cur)) / (max(cur) - min(cur))
  } else if (scale == "none") {
    # Do nothing
  } else {
    stop("Invalid scale option. Choose from 'mad', 'z', 'minmax', or 'none'.")
  }
  
  
  p0 <- ggplot(data.frame(y=cur), aes(x=y, fill = after_stat(x))) +
    theme_classic()+
    geom_histogram(bins=100) +
    xlab(paste0(gepname, " - scaled spectra score")) +
    scale_fill_viridis_c(limits=range(cur), name="GEP spectra\nscaled", option='H') +
    ggtitle(paste0(gepname, " profile")) +
    theme(plot.title = element_text(hjust = 0.5))
  plots[["hist"]] <- p0
  
  df.plot       <- data.frame(y=cur[1:k], genes=genes[1:k])
  df.plot$genes <- factor(df.plot$genes, levels=df.plot$genes)
  df.plot$gepname <- gepname
  p1 <- ggplot(df.plot, aes(x=genes, y=gepname, fill=y)) + 
    geom_tile() +
    scale_fill_viridis_c(limits=range(cur), name="GEP spectra\nscaled", option='H') +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("Top ", k," genes"))
  
  plots[["top_genes"]] <- p1
  
  if (!is.null(cyto)) {
    cur.cyto <- cur[genes %in% cyto]
    cur.cyto.gene  <- genes[genes %in% cyto]
    
    names(cur.cyto) <- cur.cyto.gene
    cur.cyto        <- sort(cur.cyto, decreasing = T)  
    
    k.cyto <- min(k, length(cur.cyto))
    
    df.plot       <- data.frame(y=cur.cyto[1:k.cyto], genes=names(cur.cyto)[1:k.cyto])
    df.plot$genes <- factor(df.plot$genes, levels=df.plot$genes)
    df.plot$gepname <- gepname
    
    p2 <- ggplot(df.plot, aes(x=genes, y=gepname, fill=y)) + 
      geom_tile() +
      scale_fill_viridis_c(limits=range(cur), name="GEP spectra\nscaled", option='H') +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0("Top ", k.cyto ," cytokine genes"))
    
    plots[["top_cyto_genes"]] <- p2
  }
  
  if (!is.null(tfs)) {
    cur.tfs <- cur[genes %in% tfs]
    cur.tfs.gene  <- genes[genes %in% tfs]
    
    names(cur.tfs) <- cur.tfs.gene
    cur.tfs        <- sort(cur.tfs, decreasing = T)  
    
    k.tfs <- min(k, length(cur.tfs))
    
    df.plot       <- data.frame(y=cur.tfs[1:k.tfs], genes=names(cur.tfs)[1:k.tfs])
    df.plot$genes <- factor(df.plot$genes, levels=df.plot$genes)
    df.plot$gepname <- gepname
    
    p3 <- ggplot(df.plot, aes(x=genes, y=gepname, fill=y)) + 
      geom_tile() +
      scale_fill_viridis_c(limits=range(cur), name="GEP spectra\nscaled", option='H') +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0("Top ", k.tfs ," TF genes"))
    
    plots[["top_TF_genes"]] <- p3
  }
  
  plot <- cowplot::plot_grid(plotlist=plots, ncol=1, rel_heights = c(2, 1, 1, 1), rel_widths=c(0.9, 1,1,1))
  return(plot)
  
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
  
  make_option(c("--threshold"), type="numeric", default=5e-5,
            help="Pvalue threshold to consider enrichments [default: %default]"),

  make_option(c("-a", "--annot"), type="character", default=NULL,
              help="Annotation file or NULL"),
              
  make_option(c("--qc"), type="character", default=NULL,
              help="QC result file or NULL"),
              
  make_option(c("--logSpectra"), action="store_true", default=FALSE,
              help="Log2-transform spectra matrix [default: %default]"),

  make_option(c("--scaleSpectra"), action="store_true", default=FALSE,
              help="Scale the spectra matrix [default: %default]"),

  make_option(c("-o", "--output"), type="character", default="output.tsv",
              help="Output prefix [default: %default]"),
              
  make_option(c("--tfFile"), type="character", default=NULL,
              help="Transcription factor file (one gene per line) [default: %default]"),

  make_option(c("--cytoFile"), type="character", default=NULL,
              help="Cytokine/marker file (one gene per line) [default: %default]")
  
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
enrich.threshold <- opt$threshold
tf.file <- opt$tfFile
cyto.file <- opt$cytoFile
qc.file <- opt$qc

scale.method <- "mad"

if (is.null(spectra.file)) {
  stop("You must provide a spectra file with -S / --spectra")
}

cat("[INFO] spectra: ", spectra.file, "\n")
cat("[INFO] enrichment: ", enrichment.file, "\n")
cat("[INFO] annotations: ", annot.file, "\n")
cat("[INFO] qc file: ", qc.file, "\n")


#-------------------------------------------------------------------------------
# Spectra
spectra     <- read(spectra.file)
spectra     <- t(spectra)

# Find the top.n genes based on the spectra file
top.spectra <- apply(apply(spectra, 2, function(x){rownames(spectra)[order(x, decreasing = T)]})[1:top.n,], 2, paste0, collapse=sep)
top.mat     <- data.frame(top_gep_genes=top.spectra)

cat("[INFO] read spectra\n")
#-------------------------------------------------------------------------------
# Create a profile overview for each gep

nrows <- 3

if (!is.null(tf.file)) {
  tfs <- fread(tf.file, data.table=F)[,1]
  nrows <-  nrows+1
} else {
  tfs <- NULL
}

if (!is.null(cyto.file)) {
  cyto <- fread(cyto.file, data.table=F)[,1]
  nrows <-  nrows+1
  
} else {
  cyto <- NULL
}

pdf(width=6, height=(nrows*1.5), file=paste0(output.prefix, "_gep_profiles.pdf"))
for (gep in 1:ncol(spectra)) {
  p <- plot.top.k.genesets(spectra, gep, cyto=cyto, tfs=tfs, scale=scale.method)
  plot(p)
}
dev.off()
#-------------------------------------------------------------------------------
# QC results

if (!is.null(qc.file)) {
  qc      <- read(qc.file)
  head(qc)
  top.mat <- cbind(top.mat, qc[rownames(top.mat),])
}

#-------------------------------------------------------------------------------
# Enrichment
if (!is.null(enrichment.file)) {
  enrich       <- read(enrichment.file, rn=F)
  cat("[INFO] read erichment\n")
  
  if (is.null(tests)) {
    tests <- unique(enrich$test)
  }
  
  if (is.null(databases)) {
    databases <- unique(enrich$database)
  }
  
  enrich       <- enrich[enrich$test %in% tests,]
  enrich       <- enrich[enrich$database %in% databases,]

  if (nrow(enrich) >=1) {
    
      # Select only the positive enrichments, not the depeletions
      enrich  <- enrich[enrich$effect_size > 0, ]
      
      enrich$group  <- paste0(enrich$test, "_", enrich$condition)
      enrich        <- enrich[,c("group", "trait", "condition", "test", "pvalue")]
      enrich        <- enrich[enrich$pvalue < enrich.threshold,]
      
      enrich$pvalue <- -log10(enrich$pvalue)
      
      # Select the top pathways
      top          <- get.topn(enrich, top.n=top.n)
      
      # Aggregate by group
      top.agg <- aggregate(top[,c("condition", "test", "trait")], by=list(top$group), function(x){paste0(unique(x), collapse = sep)})
      
      # Convert to matrix
      top.agg <- acast(top.agg, condition ~ test, value.var = "trait")
      
      # Ensure rows exactly match colnames(spectra) (GEP order).  
      # Add missing rows filled with NA and reorder as needed.
      desired_rows <- colnames(spectra)
      if (is.null(dim(top.agg))) {
        # If acast returned a vector (no columns), create empty matrix with desired rows
        top.agg <- matrix(NA, nrow=length(desired_rows), ncol=0,
                          dimnames = list(desired_rows, character(0)))
      } else {
        missing_rows <- setdiff(desired_rows, rownames(top.agg))
        if (length(missing_rows) > 0L) {
          add_mat <- matrix(NA, nrow=length(missing_rows), ncol=ncol(top.agg),
                            dimnames = list(missing_rows, colnames(top.agg)))
          top.agg <- rbind(top.agg, add_mat)
        }
        # Reorder rows to match desired order (preserves NA rows for missing entries)
        top.agg <- top.agg[desired_rows, , drop=FALSE]
      }
            
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
  
  if (scale.spectra) {
    
    scale.method="mad"
    if (scale.method=="z") {
      spectra.scaled <- scale(spectra)
    } else if (scale.method=="mad") {
      spectra.scaled <- apply(spectra, 2, function(x){(x-median(x))/mad(x)})
    } else if (scale.method=="minmax") {
      spectra.scaled <- apply(spectra, 2, function(x){(x-min(x))/(max(x)-min(x))})
    } else if (scale.method == "none") {
      # Do nothing
    } else {
      stop("Invalid scale method. Choose from 'z', 'mad', 'none' or 'minmax'.")
    }
    
    rownames(spectra.scaled) <- rownames(spectra)
    colnames(spectra.scaled) <- colnames(spectra)
  
    m               <- spectra.scaled[ol,]    
  } else {
    m               <- spectra[ol,]
  }
  
  
  #spectra.ranked <- apply(spectra, 2, function(x){rank(x, ties.method="average")})
  #spectra.ranked <- (max(spectra.ranked)+1) - spectra.ranked
  
  get_rank <- function(x, dist) {
    sum(dist >= x)
  }
  
  # Average the scores first, then find the rank of the average score.
  # This is more robust than ranking each gene and then averaging the ranks, as strong rank outliers pull down the average rank
  m2.ranked              <- aggregate(spectra[ol,], by=list(group=annot[ol,1]), FUN=mean)
  
  rownames(m2.ranked)    <- m2.ranked[,1]
  m2.ranked <- m2.ranked[,-1]
  
  for (i in 1:nrow(m2.ranked)) {
    for (j in 1:ncol(m2.ranked)) {
      m2.ranked[i,j] <- get_rank(m2.ranked[i,j], dist=spectra[,j])
    }
  }
  m2.ranked              <- t(m2.ranked)
  
  # Determine mean spectra score
  m2              <- aggregate(m, by=list(group=annot[ol,1]), FUN=mean)
  orig.grp        <- m2[,1]
  rownames(m2)    <- m2[,1]
  m2              <- t(m2[,-1])
  
  
  if (log.spectra) {
    m2 <- log2(m2+1)
  }

  m2 <- na.omit(m2)
  m2.ranked <- m2.ranked[rownames(m2),colnames(m2)]
  
  
  # Cluster 
  d <- dist(t(m2))       # transpose to cluster columns
  hc <- hclust(d)   
  
  
  orig.tm <- colnames(top.mat)

  tm.orig <- top.mat

  top.mat <- cbind(tm.orig, m2[rownames(tm.orig), hc$order])
  top.mat.ranked <- cbind(tm.orig, m2.ranked[rownames(tm.orig), hc$order])
  orig.names <- c(orig.tm, orig.grp[hc$order])  
  
  print(dim(m2))
  print(dim(m2.ranked))
  
  pdf(width=2+(nrow(m2)*1), height=2+(ncol(m2)*0.25), file=paste0(output.prefix, "_marker_heatmap.pdf"))
  plot(plot_simple_hm(m2))
  dev.off()
}

#-------------------------------------------------------------------------------
# Output
top.mat <- data.frame(rownames=rownames(top.mat), gep=paste0("GEP",rownames(top.mat)), top.mat)
colnames(top.mat) <- c("rownames", "gep", orig.names)
write.table(top.mat,
            file=paste0(output.prefix, "_summary.tsv"),
            quote=F,
            sep="\t",
            row.names=F)

  
top.mat.ranked <- data.frame(rownames=rownames(top.mat.ranked), gep=paste0("GEP",rownames(top.mat.ranked)), top.mat.ranked)
colnames(top.mat.ranked) <- c("rownames", "gep", orig.names)
write.table(top.mat.ranked,
            file=paste0(output.prefix, "_summary_ranks.tsv"),
            quote=F,
            sep="\t",
            row.names=F)
