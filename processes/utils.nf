#!/usr/bin/env nextflow

process fetch_gene_id_reference {
    label "tiny"
    scratch params.rn_scratch
    
    container params.rn_container
    conda params.rn_conda
    
    publishDir "$params.rn_publish_dir/reference/ensembl", mode: 'copy'
    
    input:
        val(version)
    output:
        path("v${version}_ensembl.tsv", emit: ensembl)
        path("v${version}_name_to_ensembl.tsv", emit: name_to_ensembl)
        path("v${version}_ensembl_to_name.tsv", emit: ensembl_to_name)
    script:
    
    cmd =
    """
    fetch_ensembl_genes.r ${version} v${version}
    """  
}

process invert_id_link {
    label "tiny"
    publishDir "$params.rn_publish_dir/reference/id_link", mode: 'copy'

    input:
        path(id_linker)
    output:
        path("*_inverted.tsv", emit: inverted)
    script:
    """
    outfile="\$(basename ${id_linker} | sed 's/.tsv//g')"
    awk '{print \$2"\t"\$1"}' ${id_linker} > \${outfile}_inverted.tsv
    """
}


process filter_biotype {
    label "tiny"
    publishDir "$params.rn_publish_dir/reference/ensembl", mode: 'copy'

    input:
        path(ensembl_file)
    output:
        path("*_biotype_filtered.ensembl.txt", emit: biotype_ensembl)
        path("*_biotype_filtered.gene_names.txt", emit: biotype_gene_name)
    script:
    
    if (params.convert.biotype_filter == null) {
        error "Biotype filter parameter is not set but process is called."
    }
    
    cmd =
    """
    # Filter ensembl file by biotype using egrep pattern from config
    cat ${ensembl_file} | \
    egrep ${params.convert.biotype_filter} | \
    awk -F'\t' '{print \$1}' > \
    \$(basename ${ensembl_file} | sed 's/.tsv/_biotype_filtered.ensembl.txt/g')
    
    # Also create a version with gene names
    cat ${ensembl_file} | \
    egrep ${params.convert.biotype_filter} | \
    awk -F'\t' '{print \$10}' > \
    \$(basename ${ensembl_file} | sed 's/.tsv/_biotype_filtered.gene_names.txt/g')
    """
    
    cmd
}


