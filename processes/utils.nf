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

process invert_id_link {
    label "tiny"
    publishDir "$params.rn_publish_dir/reference/id_link", mode: 'copy'

    input:
        path(id_linker)
    output:
        path("*_inverted.tsv", emit: inverted)
    script:
    """
    outfile="\$(basename ${id_linker} | sed 's/.tsv//g')"
    awk '{print \$2"\t"\$1"}' ${id_linker} > \${outfile}_inverted.tsv
    """
}


