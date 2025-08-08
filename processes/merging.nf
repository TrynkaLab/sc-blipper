#!/usr/bin/env nextflow

process merge_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.convert.container
    conda params.convert.conda
    
    publishDir "$params.rn_publish_dir/h5ad/merged", mode: 'symlink'
    
    input:
        val(id)
        path(files, arity: '1..*')
    output:
        tuple val(id), path("${id}_merged.h5ad")
    script:
    """
    merge_h5ad.py ${id}_merged.h5ad ${files.join(' ')}
    """
}
