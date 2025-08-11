#!/usr/bin/env nextflow

process merge_h5ad {
    label params.merge.label
    scratch params.rn_scratch
    
    container params.merge.container
    conda params.merge.conda
    
    publishDir "$params.rn_publish_dir/h5ad/merged", mode: 'symlink'
    
    input:
        val(id)
        path(files, arity: '1..*')
    output:
        tuple val(id), path("${id}_merged.h5ad"), emit: merged
        path("non_overlapping_genes.txt", optional: true, emit: non_overlaps)
    script:

    cmd = "merge_h5ad.py ${id}_merged.h5ad ${files.join(' ')}"
    
    if (params.merge.overlap) {
        cmd += " --overlap"
    }
    
    cmd
}
