#!/usr/bin/env Rscript

library(data.table)
library(readxl)



# Function to write GMT file
write_gmt <- function(gene_list, description, source, output_file) {
  write(paste(description, source, paste(gene_list, collapse="\t")), 
        file=output_file, 
        append=TRUE)
}

ensembl <- fread("ensembl/v115_ensembl_to_name.tsv", header=F, data.table = F)
rownames(ensembl) <- ensembl[,2]

url <- "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"
tmp <- tempfile(fileext = ".xls")
download.file(url, tmp, mode = "wb")


mitocarta <- read_excel(tmp, sheet="C MitoPathways")
mitocarta <- na.omit(mitocarta)[,-1]

mc.gene_name <- setNames(
  lapply(mitocarta$Genes, function(x) strsplit(x, ", ", fixed = TRUE)[[1]]),
  mitocarta$MitoPathway
)


# 7 out of 1035 genes could not be found in the gene id, remove these
length(unique(unlist(mc.gene_name))) - sum(unique(unlist(mc.gene_name)) %in% ensembl$V2)
to.keep <- unique(unlist(mc.gene_name))[unique(unlist(mc.gene_name)) %in% ensembl$V2]


mc.gene_name <- lapply(mc.gene_name, function(x){x[x %in% to.keep]})

# Convert spaces in gene names to underscores
names(mc.gene_name) <- sapply(names(mc.gene_name), function(x){gsub(" ", "_", x)})


# Convert gene name to ensembl id
mc.ensembl <- mc.gene_name
mc.ensembl <- lapply(mc.ensembl, sapply, function(x){ensembl[x,1]})


# Remove the old version
if (file.exists(paste0(file_path, "gene_sets/symbols/mitocarta.human.v3.symbols.gmt"))) {
  file.remove(paste0(file_path, "gene_sets/symbols/mitocarta.human.v3.symbols.gmt"))
}

if (file.exists(paste0(file_path, "gene_sets/ensembl/mitocarta.human.v3.ensembl.gmt"))) {
  file.remove(paste0(file_path, "gene_sets/ensembl/mitocarta.human.v3.ensembl.gmt"))
}

for (pw in names(mc.gene_name)) {
  write_gmt(mc.gene_name[[pw]], 
            pw,
            "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls",
            "gene_sets/symbols/mitocarta.human.v3.symbols.gmt")
  
}


for (pw in names(mc.ensembl)) {
  write_gmt(mc.ensembl[[pw]], 
            pw,
            "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls",
            "gene_sets/ensembl/mitocarta.human.v3.ensembl.gmt")
  
}


