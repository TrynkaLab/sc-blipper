#!/usr/bin/env nextflow

include { seurat_to_h5ad; link_h5ad; } from '../processes/convert_to_h5ad.nf'
include { merge_h5ad } from '../processes/merging.nf'
include { fetch_gene_id_reference } from "../processes/utils.nf"
include { cnmf_pre_process; cnmf_prepare; cnmf_factorize; cnmf_combine; cnmf_kselection; cnmf_consensus } from "../processes/cnmf.nf"
include { gsea; ora } from "../processes/enrichment.nf"

// Main workflow
workflow run_cnmf {
    
    main:
        // If needed, fetch ensembl file
        if (params.convert.convert_gene_names) {
            if (params.convert.id_linker != null) {
                id_linker = file(params.convert.id_linker)
                if (!id_linker.exists()) {
                    throw new Exception("Supplied id linker file does not exist")
                }
                // Custom conversion file
                id_linker = Channel.value(file(params.convert.id_linker))
            } else {                
                if (params.convert.ensembl_to_gene) {
                    // ENSEMBL > gene name conversion
                    id_linker = fetch_gene_id_reference(params.convert.ensembl_version).ensembl_to_name
                } else {
                    // gene name > ENSEMBL conversion
                    id_linker = fetch_gene_id_reference(params.convert.ensembl_version).name_to_ensembl
                }
            }
        } else {
            // In case there is no mapping file
            id_linker = Channel.value(file("NO_MAPPING"))
        }
    
        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(params.rn_manifest)
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
        
        // Converts seurat files to h5ad
        convert_out_a = seurat_to_h5ad(manifest_split.seu, id_linker).h5ad

        // Converts just symplink h5ad the files so they are in the ouput
        // In case conversion is needed, do that here
        convert_out_b = link_h5ad(manifest_split.h5ad, id_linker).h5ad

        // Merge the channels, they should now contain all h5ad files
        convert_out = convert_out_a.concat(convert_out_b)
        convert_out_flat = convert_out.flatMap{row -> {row[1]}}.collect()

        // Merge the .h5ad files into a single file
        merge_out = merge_h5ad(params.rn_runname, convert_out_flat).merged
        
        //--------------------------------------------------------
        // Optionally run cNMF preprocessing, which creates harmony corrected counts
        if (params.cnmf.preprocess) {
            // Preprocess the cnmf file
            cnmf_preprocess = cnmf_pre_process(merge_out)
            
            // Channel with run name and .h5ad file with counts
            cnmf_in = cnmf_preprocess.preproccessed
        } else {
            // Otherwise use the merged file, in case there is one batch,
            // just uses the otehr file
            cnmf_in = merge_out
        }
            
        //--------------------------------------------------------
        //                       cNMF
        //--------------------------------------------------------
        // Prepare the cnmf input folder 
        cnmf_prepared = cnmf_prepare(cnmf_in)
        
        // Create the channel with the n_workers for jobs to run in parallel
        cnmf_factorize_in = cnmf_prepared
         .map { v -> (0..(params.cnmf.n_workers-1)).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        // Perform the factorizations
        cnmf_factorize_out = cnmf_factorize(cnmf_factorize_in)
        
        // Wait for all the workers to complete, then:
        // Flatten all the output into a single path object
        cnmf_factorize_out = cnmf_factorize_out.collect()
        
        // Combine the output   
        cnmf_combine_out = cnmf_combine(cnmf_prepared, cnmf_factorize_out)     
 
        // Make the k selection plot
        cnmf_kselection(cnmf_prepared, cnmf_combine_out)
        
        // Make the consensus factorizations
        cnmf_consensus_in = cnmf_prepared
         .map { v -> (params.cnmf.k).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out)
        // This is the end of the cnmf processing

        //--------------------------------------------------------
        // Optional post-proccessing of cnmf
        
        if (params.cnmf.run_enrichment && params.enrich.gmt_files != null) {
                gmt_files = Channel.of(params.enrich.gmt_files)
            
                gsea_in = cnmf_out.spectra_k.map{i -> tuple("k_"+i[0], i[1])}
                
                gsea_in.view()
                //gsea_out = gsea("cnmf/consensus/${params.rn_runname}",
                //                gsea_in)
                
    
        }
        
        
        if (params.cnmf.run_decoupler) {
            
            
        }

}