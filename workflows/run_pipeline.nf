#!/usr/bin/env nextflow

include { seurat_to_h5ad; link_h5ad; } from '../processes/convert_to_h5ad.nf'
include { merge_h5ad } from '../processes/merging.nf'
include {fetch_gene_id_reference} from "../proccesses/convert_gene.nf"

// Main workflow
workflow run_pipeline {
    
    main:
        // If needed, fetch ensembl file
        if (params.convert.convert_gene_names) {
            id_linker = file(params.convert.ensembl_id_linker)
            
            if (!id_linker.exists()) {
                id_linker = fetch_gene_id_reference(params.convert.ensembl_version)
            }
        }
    
        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.id,
            file(row.file))}
        
        // Auto detect if it is a h5ad already or a seurat file
        manifest_split = manifest.branch{ it ->
            h5ad: it[1].extension == "h5ad"
            seu: ['rds', 'RDS', 'Rds'].contains(it[1].extension)
            other: false
        }
        
        // Converts seurat files to h5ad
        convert_out_a = seurat_to_h5ad(manifest_split.seu)

        // Converts just symplink h5ad the files so they are in the ouput
        convert_out_b = link_h5ad(manifest_split.h5ad)

        // Merge the channels, they should now contain all h5ad files
        convert_out = convert_out_a.concat(convert_out_b)
        convert_out_flat = convert_out.flatMap{row -> {row[1]}}.collect()

        // Merge the .h5ad files into a single file
        merge_out = merge_h5ad(params.rn_runname, convert_out_flat)

        // Optionally run Harmony correction as implemented by cNMF
    
}