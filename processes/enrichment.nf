#!/usr/bin/env nextlflow

process preprocess_matrix {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
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
    
    container params.rn_container
    conda params.rn_conda
    
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
        tuple val(id), path("${id}_fgsea_results.enrich.gz"), emit: std
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
    
    // Convert to standard format and zip
    cmd += 
    """
    awk -v OFS='\t' 'BEGIN {print "test","database","condition","trait","effect_size","error","effect_size_norm","pvalue"}' > ${id}_fgsea_results.enrich
    cat ${id}_fgsea_results.tsv | tail -n +2 | awk -v OFS='\t' '{print "GSEA",\$11,\$9,\$1,\$5,\$4,\$6,\$2}' >> ${id}_fgsea_results.enrich
    
    gzip -f ${id}_fgsea_results.tsv
    gzip -f ${id}_fgsea_results.enrich
    """
        
    cmd
}


process ora {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
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
        tuple val(id), path("${id}_ora_results.enrich.gz"), emit: std
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

    // Convert to standard format and zip
    cmd +=
    """
    awk -v OFS='\t' 'BEGIN {print "test","database","condition","trait","effect_size","error","effect_size_norm","pvalue"}' > ${id}_ora_results.enrich
    cat ${id}_ora_results.tsv | tail -n +2 | awk -v OFS='\t' '{print "ORA",\$10,\$8,\$1,\$4,"NA","NA",\$2}' >> ${id}_ora_results.enrich
    
    gzip ${id}_ora_results.tsv
    gzip ${id}_ora_results.enrich
    """    
    
    cmd
}


process decoupler {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/enrichments", mode: 'symlink'
    
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(file)
        val(transpose)
        path(mapping_file)
    output:
        path("${id}_*.tsv.gz", emit: decoupler)
        tuple val(id), path("${id}_*.enrich.gz"), emit: std
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
    
    // Convert to standard format and zip
    cmd +=
    """
    # Collectri
    awk -v OFS='\t' 'BEGIN {print "test","database","condition","trait","effect_size","error","effect_size_norm","pvalue"}' > ${id}_collectri.enrich
    cat ${id}_collectri_activities.tsv | tail -n +2 | awk -v OFS='\t' '{print "TF_TARGET","collecttri",\$3,\$2,\$4,"NA","NA",\$5}' >> ${id}_collectri.enrich
    
    # Progeny
    awk -v OFS='\t' 'BEGIN {print "test","database","condition","trait","effect_size","error","effect_size_norm","pvalue"}' > ${id}_progeny.enrich
    cat ${id}_progeny_activities.tsv | tail -n +2 | awk -v OFS='\t' '{print "ACTIVITY","progeny",\$3,\$2,\$4,"NA","NA",\$5}' >>  ${id}_progeny.enrich
    
    # Zip
    gzip ${id}_progeny.enrich
    gzip -f ${id}_*.tsv
    """    
    cmd
}


process concat_enrichment_results {
    label params.enrich.label
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/${prefix}/${id}/final", mode: 'symlink'
        
    input:
        // A couple of parameters are set as input values
        val(prefix)
        tuple val(id), path(files)
        path(annot_file)
    output:
        path("${id}_merged.tsv.gz", emit: std)
        path("${id}_merged_nominal.tsv", emit: std_nominal)

    script:
    
    cmd =
    """
    # Extract header from the first file and prepend 'filename' as the first column
    zcat "${files[0]}" | head -n1 | awk '{print "filename\t"\$0}' > "${id}_merged.tmp"

    for f in ${files}; do
        # Get filename without any extensions
        filename=\$(basename "\$f")
        # Strip all extensions (remove everything after the first dot)
        filename="\${filename%%.*}"
        
        # Add filename as first column, skip header line
        zcat "\$f" | tail -n +2 | awk -v fname="\$filename" '{print fname"\t"\$0}' >> "${id}_merged.tmp"
    done
    
    """
    
    if (annot_file.getFileName().toString() == "NO_ANNOT") {
        cmd += "annotate_enrichments.py --input ${id}_merged.tmp --output ${id}_merged.tsv"
    } else {
        cmd += "annotate_enrichments.py --input ${id}_merged.tmp --output ${id}_merged.tsv --annot ${annot_file}"
    }
    
    cmd +=
    """
    # Cleanup
    rm ${id}_merged.tmp
    
    colnum=$(head -1 ${id}_merged.tsv | tr '\t' '\n' | grep -n "^padj_test\$" | cut -d: -f1)
    cat ${id}_merged.tsv | awk -v c=\$colnum -F'\t' 'NR==1 || \$c < 0.05' > ${id}_merged_nominal.tsv
    
    gzip -f ${id}_merged.tsv
    """
    cmd
    
}