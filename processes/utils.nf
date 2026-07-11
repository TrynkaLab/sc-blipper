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
    awk '{print \$2"\t"\$1}' ${id_linker} > \${outfile}_inverted.tsv
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


process preprocess_matrix {
    label params.enrich.label
    scratch params.rn_scratch

    container params.rn_container
    conda params.rn_conda

    publishDir "$params.rn_publish_dir/${prefix}/${id}/processed", mode: 'symlink'

    input:
        val(prefix)
        tuple val(id), file(path), val(transpose), file(id_linker), file(universe)
        path(subset_genes)
    output:
        tuple val(id), file("${id}_processed_matrix.tsv"), emit: matrix
        path("${id}*.log", emit: log )
    script:
    cmd =
    """
    table_proccessor.py \
    --input ${path} \
    --output ${id}_processed_matrix.tmp\
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

    if (subset_genes.getFileName().toString() != "NO_SUBSET") {
        cmd +=
        """
        table_proccessor.py \
        --input ${id}_processed_matrix.tmp \
        --output ${id}_processed_matrix.tsv \
        --row-file ${subset_genes}
        """
    } else {
        cmd +=
        """
        mv ${id}_processed_matrix.tmp ${id}_processed_matrix.tsv
        """
    }

    cmd
}


// Binarize a numeric genes x conditions matrix based on top-N rank per column,
process binarize_matrix {
    label params.ldsc.label

    container params.rn_container
    conda params.rn_conda

    publishDir "$params.rn_publish_dir/${prefix}/${id}/processed", mode: 'symlink'

    input:
        val(prefix)
        val(id)
        path(matrix)
        val(top)
        val(absolute)
        val(ascending)

    output:
        path("${matrix.simpleName}__top${top}.tsv"), emit: matrix
        path("${matrix.simpleName}__top${top}.log"), emit: log

    script:
        cmd =
        """
        binarize_matrix.py \\
            --input ${matrix} \\
            --output ${matrix.simpleName}__top${top}.tsv \\
            --top ${top}\
        """

        if (absolute) {
            cmd += " --absolute"
        }

        if (ascending) {
            cmd += " --ascending"
        }

        cmd
}


// Split the genes x conditions matrix into one gene-list file per condition.
// Each output file is named {condition}.geneset.txt and contains the gene IDs
// (rows) where the column value equals 1.
process extract_conditions {
    label params.ldsc.label

    container params.rn_container
    conda params.ldsc.conda

    input:
        tuple val(tag), path(gene_matrix)

    output:
        tuple val(tag), path("*.geneset.txt"), emit: gene_sets

    script:
        cmd =
        """
        python3 << 'PYEOF'
        import pandas as pd
        import re

        matrix = pd.read_csv("${gene_matrix}", sep="\\t", index_col=0)
        for col in matrix.columns:
            genes = matrix.index[matrix[col] == 1].tolist()
            safe_col = re.sub(r'[^\\w.-]', '_', col)
            with open(f"{safe_col}.geneset.txt", "w") as fh:
                fh.write("\\n".join(genes) + "\\n")
        PYEOF
        """
        cmd
}
