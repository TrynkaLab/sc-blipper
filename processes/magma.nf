#!/usr/bin/env nextflow


process ensembl_to_magma_geneloc {
    label "tiny"
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/reference/magma", mode: 'symlink'
    
    input:
        path(ensembl_file)
        val(use_ensembl_id)
    output:
        path("*_magma_reference.gene.loc")
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    
    file='$ensembl_file'
    ensembl <- fread(file, data.table=F)
    
    ensembl[is.na(ensembl\$hgnc_symbol), "hgnc_symbol"] <- "NO_NAME"
    ensembl[ensembl\$hgnc_symbol == "",  "hgnc_symbol"] <- "NO_NAME"
    ensembl[ensembl\$hgnc_symbol == "NO_NAME",  "hgnc_symbol"] <- make.unique(ensembl[ensembl\$hgnc_symbol == "NO_NAME",  "hgnc_symbol"], sep="_")
        
    # Filter to protein coding only (like in the original column)
    ensembl <- ensembl[ensembl\$gene_biotype %in% c("protein_coding"),]
    
    # Filter to autosomal genes
    ensembl <-  ensembl[ensembl\$chromosome_name %in% c(as.character(1:22), "X", "Y"),]

    if ('$use_ensembl_id' == 'true') {
        magma <- ensembl[,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol")]
    } else {
        magma <- ensembl[,c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id")]
    }
    
    # Remove duplicated gene ids (could happen if gene names are used)
    magma <- magma[!duplicated(magma[,1]),]
    
    # Convert strand
    magma\$strand <- ifelse(magma\$strand == 1, "+", "-")
    
    out_name <- gsub(".tsv\$", "",basename(file))
    
    write.table(magma, file=paste0(out_name, "_magma_reference.gene.loc"), row.names=F, sep="\t", col.names=F, quote=F)
    """   
}


process perpare_sumstats {
    label "tiny"
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/reference/magma/sumstats", mode: 'symlink'
    
    input:
        tuple val(name), val(n), val(variant_col), val(p_col), file(file)
    output:
        tuple val(name), val(n), val(variant_col), val(p_col), path("${name}_snp_pval.tsv")
    script:
    """
    #!/usr/bin/env bash
    
    if [[ "$file" == *.gz ]]; then
        echo "File is gzipped — using zcat"
        zcat "$file" | awk -v varcol="${variant_col}" -v pcol="${p_col}" '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if (\$i == varcol) vc=i
                    if (\$i == pcol) pc=i
                }
                next
            }{print \$vc"\t"\$pc}' > ${name}_snp_pval.tsv
    
    else
        echo "File is plain text — using cat"
        cat "$file" | awk -v varcol="${variant_col}" -v pcol="${p_col}" '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if (\$i == varcol) vc=i
                    if (\$i == pcol) pc=i
                }
                next
            }{print \$vc"\t"\$pc}' > ${name}_snp_pval.tsv
    fi

    """    
}


process magma_annotate {
    label "tiny"
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/reference/magma", mode: 'symlink'
    
    input:
        val(name)
        file(snploc)
        file(geneloc)
    output:
        path("${name}.genes.annot")
    script:
    """
    magma --annotate window=${params.magma.annotate_window} \
    --snp-loc ${snploc} \
    --gene-loc ${geneloc} \
    --out ${name} \
    """
}

process magma_prepare {
    label params.magma.label
    
    container params.magma.container
    conda params.magma.conda
    
    //publishDir "$params.rn_publish_dir/magma/per_trait/${name}", mode: 'symlink'
    
    input:
        tuple val(name), val(n), val(variant_col), val(p_col), file(sumstat), val(batch)
        tuple val(plink), file(bed), file(bim), file(fam)
        file(gene_annot)
    output:
        tuple val("${name}"), path("${name}.batch*"), emit: batches
    script:
    """
    magma --bfile ${plink} \
    --pval ${sumstat} N=${n} \
    --gene-annot ${gene_annot} \
    --out ${name} \
    --batch ${batch} ${params.magma.n_batch}
    """
}


process magma_merge {
    label params.magma.label
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/magma/per_trait/${name}", mode: 'symlink'
        
    input:
        tuple val(name), path(files)
        
    output:
        tuple val("${name}"), path("${name}*.log"), emit: logs
        tuple val("${name}"), path("${name}.genes.raw"), emit: raw
        tuple val("${name}"), path("${name}.genes.out"), emit: out
        
    script:
    """
    magma --merge ${name} --out ${name}
    
    # Collect the individual logs of the batches, to diagnose any issues
    cat *batch*.log > ${name}.batches.log
    """
}


process magma_enrich {
    label params.magma.label
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/magma/per_trait/${trait}/enrich", mode: 'symlink'
        
    input:
        tuple val(trait), file(magma_raw), val(database), file(geneset)
        
    output:
        path("${database}.log", emit: log)
        path("${database}.raw", emit: raw)
        path("${database}.out", emit: out)
        
    script:
    """
    magma --gene-results ${magma_raw} \
    --set-annot ${geneset} \
    --out ${trait}_${database} \
    """
}


process magma_assoc {
    label params.magma.label
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/magma/per_trait/${trait}/assoc", mode: 'symlink'
        
    input:
        tuple val(trait), file(magma_raw), val(database), file(matrix)
        val(transpose)
        
    output:
        path("${trait}_${database}.gsa.out", emit: out)
        path("${trait}_${database}.log", emit: logs)
    script:
    
    cmd =
    """
    if [[ "$matrix" == *.gz ]]; then
        echo "File is gzipped, unzipping"
        gunzip -c $matrix > unzipped.tsv
    
    else
        echo "File is already unzipped, linking"
        ln -s $matrix unzipped.tsv
    fi
    
    input="unzipped.tsv"
    """
    
    if (transpose) {
    cmd += 
    """
    > transposed.tsv 
    n_cols=\$(head -2 "unzipped.tsv" | tail -n 1 | awk -F'\t' '{print NF}')
    for i in \$(seq 1 \$n_cols); do
            
        # If it is the first column, and the rownames id is missing, add it here
        if (( i == 1 )); then
            cut -f "\$i" "unzipped.tsv" | sed "s/^\$/ROW/g" | paste -s -d '\t' >> transposed.tsv
        else
            cut -f "\$i" "unzipped.tsv" | paste -s -d '\t' >> transposed.tsv
        fi
    done
    rm unzipped.tsv
    input="transposed.tsv"
    """
    }
    
    cmd +=
    """
    magma --gene-results ${magma_raw} \
    --gene-covar \$input \
    --out ${trait}_${database} \
    """
    
    cmd
}



process magma_concat {
    label params.magma.label
    
    container params.magma.container
    conda params.magma.conda
    
    publishDir "$params.rn_publish_dir/magma/", mode: 'symlink'
        
    input:
        val(name)
        path(files)
    output:
        path("${name}_merged_magma_results.tsv", emit: out)
        
    script:
    """
    concat_magma_output.py ${name}_merged_magma_results.tsv ${files}
    """
}