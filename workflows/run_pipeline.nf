#!/usr/bin/env nextflow

include { seurat_to_h5ad; link_h5ad; } from '../processes/convert_to_h5ad.nf'

// Main workflow
workflow run_pipeline {
    
    main:
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
        convert_out = convert_out_a.merge(convert_out_b)
        
        //convert_out.view()
    
        //module load HGI/softpack/users/cc53/R_base/1
    
}