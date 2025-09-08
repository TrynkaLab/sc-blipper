#!/usr/bin/env nextflow


process seurat_to_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/h5ad/per_batch/${params.rn_runname}", mode: 'symlink'
    
    input:
        tuple val(id), path(file), val(convert_ids)
        path(mapping_file)
        path(subset_file)
    output:
        tuple val(id), path("${id}.h5ad"), emit: h5ad
        path("${id}_unmapped_genes.txt", optional: true)
    script:
    
        if (mapping_file.getFileName().toString() == "NO_MAPPING" || !convert_ids) {
            cmd = "seurat_to_h5ad.r -i ${file} -o ${id}.h5ad"
        } else {
            cmd = "seurat_to_h5ad.r -i ${file} -o ${id}.h5ad -m ${mapping_file} -u ${id}_unmapped_genes.txt"
        }
        
        if (subset_file.getFileName().toString() != "NO_SUBSET") {
            cmd += " -f ${subset_file}"
        }
        cmd
}


process link_h5ad {
    label params.convert.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/h5ad/per_batch/${params.rn_runname}", mode: 'symlink'
    
    input:
        tuple val(id), path(file, name: "input.h5ad"), val(convert_ids)
        path(mapping_file)
        path(subset_file)
    output:
        tuple val(id), path("${id}.h5ad"), emit: h5ad
        path("${id}_unmapped_genes.txt", optional: true)
    script: 
        if ((mapping_file.getFileName().toString() == "NO_MAPPING" || !convert_ids) && subset_file.getFileName().toString() == "NO_SUBSET") {
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
            
            if ((mapping_file.getFileName().toString() == "NO_MAPPING" || !convert_ids)) {
                cmd =
                """
                echo '[INFO] Subsetting genes for h5'

                update_gene_ids.py -i ${file} -o ${id}.h5ad -f ${subset_file}
                """
            } else {
                cmd  = "echo '[INFO] Updating gene ids for h5'\n\n"
                
                cmd += "update_gene_ids.py -i ${file} -m ${mapping_file} -u ${id}_unmapped_genes.txt -o ${id}.h5ad"
                
                if (subset_file.getFileName().toString() != "NO_SUBSET") {
                    cmd += " -f ${subset_file}"
                }
            }
        }
        
        cmd
    
}