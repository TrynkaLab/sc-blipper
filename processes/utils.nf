#!/usr/bin/env nextflow

process fetch_gene_id_reference {
    label "tiny"
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/reference/ensembl", mode: 'copy'
    
    input:
        val(version)
    output:
        path("v${version}_ensembl.tsv", emit: ensembl)
        path("v${version}_name_to_ensembl.tsv", emit: name_to_ensembl)
        path("v${version}_ensembl_to_name.tsv", emit: ensembl_to_name)
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


