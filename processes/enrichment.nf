#!/usr/bin/env nextlflow



process gsea {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}", mode: 'symlink'
    
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(file)
        path(gmt_files)
        path(universe)
        val(transpose)
    output:
        path("${id}_fgsea_results.tsv.gz")
    script:
    
    cmd =
    """
    run_gsea.r \
    --matrix ${file} \
    --gmt ${gmt_files.join(',')}
    """
    
    if (transpose) {
        cmd += " --transpose"
    }
    
    if (universe.getFileName().toString() != "NO_UNIVERSE") {
        cmd += " --universe ${universe}"
    }
    
    cmd += "gzip ${id}_fgsea_results.tsv"
    
    cmd
}


process ora {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.enrich.container
    conda params.enrich.conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}", mode: 'symlink'
    
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
        path("${id}_ora_results.tsv.gz")
    script:
    
    cmd =
    """
    run_gsea_ora.r \
    --matrix ${file} \
    --gmt ${gmt_files.join(',')} \
    --threshold ${threshold}
    """
    
    if (transpose) {
        cmd += " --transpose"
    }
    
    if (universe.getFileName().toString() != "NO_UNIVERSE") {
        cmd += " --universe ${universe}"
    }
    
    if (absolute) {
        cmd += " --absolute"
    }
    
    if (use_top != null) {
        cmd += " --use_top ${use_top.join(',')}"
    }
    
    cmd += "gzip ${id}_ora_results.tsv"
    
    cmd
}