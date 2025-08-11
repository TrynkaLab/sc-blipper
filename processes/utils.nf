#!/usr/bin/env nextflow

process fetch_gene_id_reference {
    label "tiny"
    scratch params.rn_scratch
    
    container params.convert.container
    conda params.convert.conda
    
    publishDir "$params.rn_publish_dir/reference/ensembl", mode: 'copy'
    
    input:
        val(version)
    output:
        path("v*_ensembl.tsv", emit: ensembl)
        path("v*_name_to_ensembl.tsv", emit: name_to_ensembl)
        path("v*_ensembl_to_name.tsv", emit: ensembl_to_name)
    script:
    
    cmd =
    """
    fetch_ensembl_genes.r ${version} v${version}
    """  
}


