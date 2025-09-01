#!/usr/bin/env nextflow


process ensembl_to_magma_geneloc {
    label "tiny"
    
    container params.rn_container
    conda params.rn_conda
    
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

    if ('$params.magma.remove_hla_genes' == 'true') {
        # Remove genes on chr6:25,726,063-33,400,644 (GRCh38)
        ensembl <- ensembl[ !(ensembl\$chromosome_name == "6" &
                                ensembl\$start_position <= 33400644 &
                                ensembl\$end_position >= 25726063), ]
    }

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

    if ('$params.magma.remove_hla_genes' == 'true') {
        write.table(magma, file=paste0(out_name, "_hla_rm_magma_reference.gene.loc"), row.names=F, sep="\t", col.names=F, quote=F)
    } else {
        write.table(magma, file=paste0(out_name, "_magma_reference.gene.loc"), row.names=F, sep="\t", col.names=F, quote=F)
    }
    """   
}


process perpare_sumstats {
    label "tiny"
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/reference/magma/sumstats", mode: 'symlink'
    
    input:
        tuple val(name), val(n), val(variant_col), val(p_col), file(file)
    output:
        tuple val(name), val(n), val(variant_col), val(p_col), path("${name}_snp_pval.tsv")
    script:
    """
    parse_sumstats.py --input ${file} --snp-col ${variant_col} --pval-col ${p_col} --output ${name}_snp_pval.tsv
    """    
}

process magma_annotate {
    label "tiny"
    
    container params.rn_container
    conda params.rn_conda
    
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
    
    container params.rn_container
    conda params.rn_conda
    
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
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/reference/magma/gene_scores/${name}", mode: 'symlink'
        
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
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/magma/${params.rn_runname}/${trait}/enrich", mode: 'symlink'
        
    input:
        tuple val(trait), file(magma_raw), val(database), file(geneset)
        path(universe)
        path(id_linker)
    output:
        path("${trait}__${database}.gsa.out", emit: out)
        path("${trait}__${database}.log", emit: logs)
        tuple val("${database}"), path("${trait}__${database}.gsa.out"), emit: per_database
        
    script:
    cmd=
    """
    magma --gene-results ${magma_raw} \
    --set-annot ${geneset} \
    --out ${trait}__${database} \
    --settings abbreviate=0\
    """
    
    if (universe.getFileName().toString() != "NO_UNIVERSE") { 
        cmd += " gene-include=${universe}"
    }
    
}

process magma_assoc {
    label params.magma.label
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/magma/${params.rn_runname}/${trait}/assoc", mode: 'symlink'
        
    input:
        tuple val(trait), file(magma_raw), val(database), file(matrix)
        val(transpose)
        path(universe)
        path(id_linker)
    output:
        path("${trait}__${database}.gsa.out", emit: out)
        path("${trait}__${database}.log", emit: logs)
        tuple val("${database}"), path("${trait}__${database}.gsa.out"), emit: per_database
    script:
    cmd =
    """
    table_proccessor.py \
    --input ${matrix} \
    --output magma_prepared.tsv\
    """
    
    if (transpose) { 
        cmd += " --transpose"
        if (universe.getFileName().toString() != "NO_UNIVERSE") { 
            cmd += " --col-file ${universe}"
        }  
        if (id_linker.getFileName().toString() != "NO_MAPPING") { 
            cmd += " --update-cols ${id_linker}"
        }
    } else {
        if (universe.getFileName().toString() != "NO_UNIVERSE") { 
            cmd += " --row-file ${universe}"
        }  
        if (id_linker.getFileName().toString() != "NO_MAPPING") { 
            cmd += " --update-rows ${id_linker}"
        }
    }
    
    cmd +=
    """
    magma --gene-results ${magma_raw} \
    --gene-covar magma_prepared.tsv \
    --settings abbreviate=0 \
    --out ${trait}__${database}
    
    # Cleanup
    rm magma_prepared.tsv
    """
    
    cmd
}

process magma_concat {
    label params.magma.label
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/enrichments", mode: 'symlink'
        
    input:
        val(prefix)
        tuple val(id), path(files)
    output:
        path("${id}_merged_magma_results.tsv.gz", emit: out)
        tuple val(id), path("${id}_merged_magma_results.enrich.gz"), emit: std

    script:
    """
    concat_magma_output.py ${id}_merged_magma_results.tsv ${files}
    
    # Convert to standard format and zip
    awk -v OFS='\t' 'BEGIN {print "test","database","condition","trait","effect_size","error","effect_size_norm","pvalue"}' > ${id}_merged_magma_results.enrich
    cat ${id}_merged_magma_results.tsv | tail -n +2 | awk -v OFS='\t' '{print "MAGMA",\$10,\$1,\$9,\$4,\$6,\$5,\$7}' >> ${id}_merged_magma_results.enrich
    
    # Zip
    gzip ${id}_merged_magma_results.tsv
    gzip ${id}_merged_magma_results.enrich
    """
}

process magma_write_manifest {
    label "tiny"
    
    publishDir "$params.rn_publish_dir/reference/magma/", mode: 'symlink'
    
    input:
        val(rows)
    output:
        file("*_magma_manifest.tsv")
    script:
    """
    echo -e "name\traw" > ${params.rn_runname}_magma_manifest.tsv
    echo -e "${rows.join('\\n')}" >> ${params.rn_runname}_magma_manifest.tsv
    """
}