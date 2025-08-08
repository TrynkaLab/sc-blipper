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
        path("linker_*")
    script:
    
    //if (file.fileName.name == "NO_LINKER") {
    cmd =
    """
    fetch_ensembl_genes.r ${version} linker_ensembl_v${version}.tsv
    """  
    //} else {
    //    cmd =
    //    """
    //    mv ${file} linker_${file.fileName.name}
    //    """
    // }
}

