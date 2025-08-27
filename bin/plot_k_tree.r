#!/usr/bin/env Rscript

library(data.table)
library(igraph)
library(ggraph)
library(optparse)

#-------------------------------------------------------------------------------
# Function to read cnmf files
read.cnmf <- function(file) {
  data <- fread(file, data.table=F)
  rownames(data) <- data[,1]
  data <- data[,-1]
  data <- t(data)
  return(data)
}

# Convert long to matrix format
to.mat <- function(flat_mat) {
  # Determine dimensions
  n_rows <- length(unique(flat_mat[,1]))
  n_cols <- length(unique(flat_mat[,2]))
  
  # Initialize empty matrix
  mat <- matrix(NA, nrow = n_rows, ncol = n_cols)
  rownames(mat) <- unique(flat_mat[,1])
  colnames(mat) <- unique(flat_mat[,2])
  
  # Fill the matrix using indexing
  mat[cbind(flat_mat[,1], flat_mat[,2])] <- flat_mat[,3]
  
  return(mat)
}

# Make the grah plot
make.plot <- function(edges, mode, order=F) {
  # Assuming 'edges' is your edge dataframe and 'g' is the graph object
  edges <- res.cur[,c("ref_gep_id", "qry_gep_id", "beta")]
  
  # Combine and get unique nodes overall
  nodes <- unique(c(unique(edges$ref_gep_id), unique(edges$qry_gep_id)))
  
  # Convert to data frame with a 'name' column expected by igraph
  nodes <- data.frame(name = nodes, group=gsub("_\\d+", "", nodes), stringsAsFactors = FALSE)
  
  # Create the edgelist
  g <- graph_from_data_frame(d = na.omit(edges), vertices = nodes, directed = T)
  
  # Create initial ggraph layout with tree layout
  p <- ggraph(g, layout = "tree") 
  
  # Extract layout data
  layout.df <- data.frame(x=NA, y=NA, group=nodes$group)
  label.df  <- data.frame()
  
  if (order) {
    layout.df <- data.frame(p$data)
    y <- 0
    for (cur.group in unique(nodes$group)) {
      selector <- layout.df$group == cur.group
      
      rk <- rank(layout.df[selector, "x"], na.last = T, ties.method="first")
      layout.df[selector, "y"] <- y 
      layout.df[selector, "x"] <- rk - mean(rk)
      label.df <- rbind(label.df, c(0, y, cur.group))
      y <- y + 1
    }

  } else {
    y <- 0
    for (cur.group in unique(nodes$group)) {
      selector <- layout.df$group == cur.group
      layout.df[selector, "y"] <- y 
      layout.df[selector, "x"] <- ((1:sum(selector)) - mean(1:sum(selector)))
      label.df <- rbind(label.df, c(0, y, cur.group))
      y <- y + 1
    }
  }

  layout.df$y <- max(layout.df$y)-layout.df$y 
  
  colnames(label.df) <- c("x", "y", "label")
  label.df$x <- min(layout.df$x) -3
  label.df$y <- as.numeric(label.df$y)
  label.df$y <- max(label.df$y) - label.df$y
  
  if (mode == "lm") {
    colscale.name <- "beta"
    widthscale.name <- "abs(beta)"
  } else if (mode == "cor") {
    colscale.name <- "Pearson R"
    widthscale.name <- "|Pearson R|"
  }
  
  # Plot with the adjusted layout
  p1 <- ggraph(g, layout = "manual", x = layout.df$x, y = layout.df$y) +
    geom_edge_link(aes(edge_alpha = abs(beta),
                       edge_width = abs(beta),
                       edge_colour = beta),
                   show.legend = TRUE) +
    geom_node_point(size = 10, shape=21, fill = "white", stroke=1, color="black" ) +
    geom_node_text(aes(label = gsub("k\\d+_", "", name)), vjust = 0.5, repel = FALSE) +
    scale_edge_color_gradient2(low = "#440154FF", mid = "grey90", high = "#F50029", midpoint = 0, name = colscale.name) +
    scale_edge_width(name = widthscale.name) +
    scale_edge_alpha(name = widthscale.name) +
    geom_text(aes(x=x, y=y, label=label), data=label.df) +
    theme_void()
  
  nnode.max <- max(aggregate(nodes[,2], by=list(nodes$group), length)[,2])
  
  return(list(plot=p1, nodes=nodes, edges=edges, ngroup=length(unique(nodes$group)), nnode.max=nnode.max))
}

#-------------------------------------------------------------------------------
#                                  Argparse
#-------------------------------------------------------------------------------
# Define options
option_list <- list(
  make_option(c("-q", "--qry"), type="character", default=NULL,
              help="Query as '1=file1 2=file2 ...'"),
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="Reference file (optional)"),
  make_option(c("-t", "--threshold"), type="character", default='auto',
              help="Threshold: numeric value, 'auto', 'signif', or NULL [default %default]"),
  make_option(c("-m", "--mode"), type="character", default="cor",
              help="Mode: 'lm' or 'cor' [default %default]"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="Output prefix [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# threshold
threshold <- opt$threshold
if (!is.null(threshold)) {
  if (threshold %in% c("auto", "signif", "NULL")) {
    if (threshold == "NULL") {
      threshold <- NULL
    }
    # leave as string for "auto" or "signif"
  } else if (suppressWarnings(!is.na(as.numeric(threshold)))) {
    threshold <- as.numeric(threshold)
  } else {
    stop("Invalid value for --threshold: must be numeric, 'auto', 'signif', or NULL")
  }
}

# mode
mode <- match.arg(opt$mode, choices=c("lm", "cor"))

# ref
ref <- opt$ref

# output
output <- opt$output

# Handle qry → Parse into named list
qry <- NULL
if (!is.null(opt$qry)) {
  # Split by space
  parts <- unlist(strsplit(opt$qry, "\\s+"))
  key_val <- strsplit(parts, "=")
  keys <- vapply(key_val, function(x) x[1], character(1))
  vals <- vapply(key_val, function(x) x[2], character(1))
  qry <- as.list(vals)
  names(qry) <- keys
}

# Print results for testing
cat("Parsed options:\n")
cat("Threshold:", threshold, "\n")
cat("Mode:", mode, "\n")
cat("Reference:", opt$ref, "\n")
cat("Query list:\n")
print(qry)

#-------------------------------------------------------------------------------
#                              Main code loop
#-------------------------------------------------------------------------------
# Sort by increasing k, this is important as the logic relies on this!
qry <- qry[as.character(sort(as.numeric(names(qry))))]

if (is.null(ref)) {
  ref.cnmf.orig     <- read.cnmf(qry[[1]])
  ref.gep.name.orig <- paste0("k",names(qry)[[1]])
  qry               <- qry[-1]
  iterative         <- TRUE
  
} else {
  ref.cnmf.orig     <- read.cnmf(ref)
  ref.gep.name.orig <- "ref"
  iterative         <- FALSE
}

cn           <- c("file", "group", "ref_gep","ref_gep_id", "qry_gep","qry_gep_id", "beta", "se", "tstat", "pval")
res          <- NULL
ref.cnmf     <- ref.cnmf.orig
ref.gep.name <- ref.gep.name.orig

for (qry.name in names(qry)) {
  cat("[INFO] Running ", qry.name, "\n")
  qry.cnmf <- read.cnmf(qry[[qry.name]])
  ol       <- intersect(rownames(ref.cnmf), rownames(qry.cnmf))

  for (i in 1:ncol(ref.cnmf)) {
    if (mode == "lm") {
      # Calculate the coefficients of a multiple regression model
      m <- lm(ref.cnmf[ol, i] ~., data=data.frame(qry.cnmf[ol,]))
      
      c           <- summary(m)$coefficients[-1,]
      rownames(c) <- gsub("X", "", rownames(c))
      
      c.f <- cbind(data.frame(file=basename(qry[[qry.name]]),
                              group=qry.name,
                              ref_gep=i,
                              ref_gep_id=paste0(ref.gep.name, "_", i),
                              qry_gep=rownames(c),
                              qry_gep_id=paste0("k",qry.name, "_",rownames(c))),
                   c)
      colnames(c.f) <- cn
      
      if (is.null(res)) {
        res <- c.f
      } else {
        res <- rbind(res, c.f)
      }
      
    } else if (mode == "cor") {
      # Calculate the individual pearson correlations instead
      for (j in 1:ncol(qry.cnmf)) {
        m           <- lm(scale(ref.cnmf[ol, i]) ~ -1+ ., data=data.frame(scale(qry.cnmf[ol,j, drop=F])))
        c           <- summary(m)$coefficients
        rownames(c) <- gsub("X", "", rownames(c))
        
        c.f <- cbind(data.frame(file=basename(qry[[qry.name]]),
                                group=qry.name,
                                ref_gep=i,
                                ref_gep_id=paste0(ref.gep.name, "_", i),
                                qry_gep=j,
                                qry_gep_id=paste0("k",qry.name, "_",rownames(c))),
                     c)
        colnames(c.f) <- cn
        
        if (is.null(res)) {
          res <- c.f
        } else {
          res <- rbind(res, c.f)
        }
      }
    }
  }
  
  if (iterative) {
    ref.cnmf     <- qry.cnmf
    ref.gep.name <- paste0("k",qry.name)
  }
}

rm(qry.cnmf, ref.cnmf, m)
gc()

#-------------------------------------------------------------------------------
if (iterative) {
  res.tmp <- res
  res.tmp$group <- 1
} else {
  res.tmp <- res
} 


for (group in unique(res.tmp$group)) {
  
  res.cur <- res.tmp[res.tmp$group == group,]
  
  # Save the edgelist
  write.table(res.cur, sep="\t", row.names=F, quote=F, file=paste0(output, "_edgelist.tsv"))
  
  # Filter edges
  if (threshold == "signif") {
    res.cur[p.adjust(res.cur$pval, method="bonf") > 0.05,"beta"] <- NA
  } else if (threshold == "auto") {
    res.cur[abs(res.cur$beta) < quantile(abs(res.cur$beta), probs=0.9),"beta"] <- NA
  } else if (is.numeric(threshold)) {
    res.cur[abs(res.cur$beta) < threshold,"beta"] <- NA
  } else {
    stop("Provided value for threshold is not valid")
  }
  
  pres <- make.plot(res.cur, mode=mode, order=F)
  pdf(width=6+(pres$nnode.max*0.25), height=5+(pres$ngroup*1), file=paste0(output, "_", group, "_ordered.pdf"))
  plot(pres$plot)
  dev.off()
  
  pres <- make.plot(res.cur, mode=mode, order=T)
  pdf(width=6+(pres$nnode.max*0.25), height=5+(pres$ngroup*1), file=paste0(output, "_", group, "_layout.pdf"))
  plot(pres$plot)
  dev.off()
}




