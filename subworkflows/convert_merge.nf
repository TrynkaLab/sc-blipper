#!/usr/bin/env nextflow
include { seurat_to_h5ad; link_h5ad; } from '../processes/convert_to_h5ad.nf'
include { merge_h5ad } from '../processes/merging.nf'

workflow convert_and_merge {

    take:
        rn_manifest
        rn_runname
        id_linker
    main:

        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.id,
            file(row.file),
            row.convert_ids.equalsIgnoreCase("true"))}
        
        // Auto detect if it is a h5ad already or a seurat file
        manifest_split = manifest.branch{ it ->
            h5ad: it[1].extension == "h5ad"
            seu: ['rds', 'RDS', 'Rds'].contains(it[1].extension)
            other: false
        }
        
        if (params.convert.subset_genes != null){
            subset_file=Channel.value(file(params.convert.subset_genes))
        } else {
            subset_file=Channel.value(file("NO_SUBSET"))
        }
        
        // Converts seurat files to h5ad
        convert_out_a = seurat_to_h5ad(manifest_split.seu, id_linker, subset_file).h5ad

        // Converts just symplink h5ad the files so they are in the ouput
        // In case conversion is needed, do that here
        convert_out_b = link_h5ad(manifest_split.h5ad, id_linker).h5ad

        // Merge the channels, they should now contain all h5ad files
        convert_out = convert_out_a.concat(convert_out_b)
        convert_out_flat = convert_out.flatMap{row -> {row[1]}}.collect()

        // Merge the .h5ad files into a single file
        merge_out = merge_h5ad(rn_runname, convert_out_flat).merged
    emit:
        merge_out = merge_out

}