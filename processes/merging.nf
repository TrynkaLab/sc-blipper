#!/usr/bin/env nextflow

process merge_h5ad {
    label params.merge.label
    scratch params.rn_scratch

    container params.rn_container
    conda params.rn_conda

    publishDir "$params.rn_publish_dir/h5ad/merged/${params.rn_runname}", mode: 'symlink'

    input:
        val(id)
        path(files, arity: '1..*')

    output:
        tuple val(id), path("${id}_merged.h5ad"), emit: merged
        path("non_overlapping_genes.tsv", optional: true, emit: non_overlaps)

    script:
    
        if (files.size() == 1) {
            cmd = 
            """
            # Only one file copy it to the expected output name
            cp ${files[0]} ${id}_merged.h5ad
            """
        } else {
            cmd = "merge_h5ad.py ${id}_merged.h5ad ${files.join(' ')}"
            if (params.merge.overlap_genes) {
                cmd += " --overlap"
            }
        }
        
        cmd
}
