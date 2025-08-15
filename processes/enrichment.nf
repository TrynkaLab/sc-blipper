#!/usr/bin/env nextlflow

process preprocess_matrix {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/processed", mode: 'symlink'

    input:
        val(prefix)
        tuple val(id), file(path), val(transpose), file(id_linker), file(universe)
    output:
        tuple val(id), file("${id}_processed_matrix.tsv"), emit: matrix
        path("${id}_processed_matrix.log", emit: log )
    script:
    cmd =
    """
    table_proccessor.py \
    --input ${path} \
    --output ${id}_processed_matrix.tsv\
    """
    
    if (transpose) { 
        cmd += " --transpose"
        if (universe.getFileName().toString() != "NO_UNIVERSE") { 
            cmd += " --col-file ${universe}"
        }  
        if (id_linker.getFileName().toString() != "NO_MAPPING") { 
            cmd += " --update-cols ${id_linker}"
        }
    } else {
        if (universe.getFileName().toString() != "NO_UNIVERSE") { 
            cmd += " --row-file ${universe}"
        }  
        if (id_linker.getFileName().toString() != "NO_MAPPING") { 
            cmd += " --update-rows ${id_linker}"
        }
    }    
    
    cmd
}


process gsea {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/enrichments", mode: 'symlink'
    
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(file)
        path(gmt_files)
        path(universe)
        val(transpose)
    output:
        path("${id}_fgsea_results.tsv.gz")
        path("${id}_fgsea_results_top*.tsv", emit: top)
        path("*.pdf", emit: plots)
    script:
    
    cmd =
    """
    run_gsea.r \
    --output_prefix ${id} \
    --matrix ${file} \
    --gmt ${gmt_files.join(',')} \
    """
    
    if (transpose) {
        cmd += " --transpose TRUE"
    }
    
    if (universe.getFileName().toString() != "NO_UNIVERSE") {
        cmd += " --universe ${universe}"
    }
    
    cmd += "\n\ngzip -f ${id}_fgsea_results.tsv"
    
    cmd
}


process ora {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/enrichments", mode: 'symlink'
    
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(file)
        path(gmt_files)
        path(universe)
        val(transpose)
        val(absolute)
        val(threshold)
        val(use_top)
    output:
        path("${id}_ora_results.tsv.gz", emit: all)
        path("${id}_ora_results_top*.tsv", emit: top)
        path("*.pdf", emit: plots)
    script:
    
    cmd =
    """
    run_gsea_ora.r \
    --output_prefix ${id} \
    --matrix ${file} \
    --gmt ${gmt_files.join(',')} \
    --threshold ${threshold} \
    """
    
    if (transpose) {
        cmd += " --transpose TRUE"
    }
    
    if (universe.getFileName().toString() != "NO_UNIVERSE") {
        cmd += " --universe ${universe}"
    }
    
    if (absolute) {
        cmd += " --absolute TRUE"
    }
    
    if (use_top != null) {
        cmd += " --use_top ${use_top.join(',')}"
    }
    
    cmd += "\n\ngzip -f ${id}_ora_results.tsv"
    
    cmd
}


process decoupler {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/enrichments", mode: 'symlink'
    
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(file)
        val(transpose)
        path(mapping_file)
    output:
        path("${id}_*.tsv.gz", emit: decoupler)
        path("*.pdf", emit: plots)
    script:
    
    cmd =
    """
    run_decoupler.r \
    --output_prefix ${id} \
    --matrix ${file}\
    """
    
    if (transpose) {
        cmd += " --transpose TRUE"
    }
    
    if (mapping_file.getFileName().toString() != "NO_MAPPING") {
        cmd += " --id_linker ${mapping_file}"
    }
    
    if (params.enrich.omnipath_cache_dir != null) {
        cmd += " --cache_dir ${params.enrich.omnipath_cache_dir}"
    }    
    
    cmd += "\n\ngzip -f ${id}_*.tsv"
    
    cmd
}