#!/usr/bin/env nextflow

process seurat_to_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.convert.container
    conda params.convert.conda
    
    publishDir "$params.rn_publish_dir/h5ad", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
    output:
        tuple val(id), path("${id}.h5ad")
    script:
    """
    seurat_to_h5ad.r ${file} ${id}.h5ad
    """
}


process link_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    publishDir "$params.rn_publish_dir/h5ad", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
    output:
        tuple val(id), path("${id}.h5ad")
    script:
    """
    if [ "${file}" != "${id}.h5ad" ]; then
        mv "${file}" "${id}.h5ad"
    else
        echo "Source and destination are the same. Skipping move."
    fi
    """
    
}