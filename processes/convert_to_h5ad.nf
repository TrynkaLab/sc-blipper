#!/usr/bin/env nextflow


process seurat_to_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/h5ad/per_batch", mode: 'symlink'
    
    input:
        tuple val(id), path(file), val(convert_ids)
        path(mapping_file)
    output:
        tuple val(id), path("${id}.h5ad"), emit: h5ad
        path("${id}_unmapped_genes.txt", optional: true)
    script:
    
        if (mapping_file.getFileName().toString() == "NO_MAPPING" || !convert_ids) {
            cmd = "seurat_to_h5ad.r ${file} ${id}.h5ad"
        } else {
            cmd = "seurat_to_h5ad.r ${file} ${id}.h5ad ${mapping_file} ${id}_unmapped_genes.txt"
        }
        cmd
}


process link_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/h5ad/per_batch", mode: 'symlink'
    
    input:
        tuple val(id), path(file, name: "input.h5ad"), val(convert_ids)
        path(mapping_file)
    output:
        tuple val(id), path("${id}.h5ad"), emit: h5ad
        path("${id}_unmapped_genes.txt", optional: true)
    script: 
        if (mapping_file.getFileName().toString() == "NO_MAPPING" || !convert_ids) {
            cmd  =
            """
            echo "[INFO] Nothing to do, just renaming the link"
            
            if [ "${file}" != "${id}.h5ad" ]; then
                cp "${file}" "${id}.h5ad"
            else
                echo "Source and destination are the same. Skipping move."
            fi
            """
        } else {
            cmd  =
            """
            echo "[INFO] Updating gene ids for h5"
            
            update_gene_ids.py ${file} ${mapping_file} ${id}_unmapped_genes.txt ${id}.h5ad
            """ 
        }
        
        cmd
    
}