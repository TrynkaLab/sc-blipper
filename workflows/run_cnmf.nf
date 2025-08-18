#!/usr/bin/env nextflow

// Processes
include { cnmf_pre_process; cnmf_prepare; cnmf_factorize; cnmf_combine; cnmf_kselection; cnmf_consensus } from "../processes/cnmf.nf"
include { gsea; ora; decoupler } from "../processes/enrichment.nf"
include { magma_assoc } from "../processes/magma.nf"
include { magma_concat as magma_concat_main } from "../processes/magma.nf"
include { magma_concat as magma_concat_per_k } from "../processes/magma.nf"
include { concat_enrichment_results as concat_enrichment_results_main } from "../processes/enrichment.nf"
include { concat_enrichment_results as concat_enrichment_results_per_k } from "../processes/enrichment.nf"

// Subworkflows
include { fetch_id_linker } from "../subworkflows/id_linking.nf"
include { convert_and_merge } from "../subworkflows/convert_merge.nf"
include { magma_base } from "../subworkflows/magma.nf"

// Main workflow
workflow cnmf {
    //TODO: make gene linker file ids unique
    main:

        //------------------------------------------------------------
        // Id linking and ensembl reference
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        
        // Add the output in the global space for convenience
        id_linker         = fetch_id_linker.out.id_linker
        id_linker_inv     = fetch_id_linker.out.id_linker_inv
        ensembl_reference = fetch_id_linker.out.ensembl_reference

        if (params.convert.convert_gene_names) {
            converter = id_linker
        } else {
            converter = Channel.value("NO_MAPPING")
        }   
        //------------------------------------------------------------ 
        // Merge seurat and H5 files
        convert_and_merge(params.rn_manifest, params.rn_runname, converter)
        
        merge_out = convert_and_merge.out.merge_out
        
        //--------------------------------------------------------
        // Optionally run cNMF preprocessing, which creates harmony corrected counts
        if (params.cnmf.preprocess) {
            // Preprocess the cnmf file
            cnmf_preprocess = cnmf_pre_process(merge_out)
            
            // Channel with run name and .h5ad file with counts
            cnmf_in = cnmf_preprocess.preproccessed
            cnmf_in_tpm = cnmf_preprocess.normalized
            cnmf_in_hvg = cnmf_preprocess.hvg
        } else {
            // Otherwise use the merged file, in case there is one batch,
            // just uses the other file
            cnmf_in = merge_out
            cnmf_in_tpm = Channel.value(file("NO_TPM"))
            cnmf_in_hvg = Channel.value(file("NO_HVG"))
        }
            
        //--------------------------------------------------------
        //                       cNMF
        //--------------------------------------------------------
        // Prepare the cnmf input folder 
        cnmf_prepared = cnmf_prepare(cnmf_in, cnmf_in_tpm, cnmf_in_hvg)
        
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
         .map { v -> (params.cnmf.k.split(",")).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out)
        // This is the end of the cnmf processing

        //--------------------------------------------------------
        // Optional post-proccessing of cnmf
        //--------------------------------------------------------

        //----------------------------------------------------------------------------------
        // Run gene set enrichment
        if (params.cnmf.run_enrichment) {
            
            // Channel with gmt files
            if (params.enrich.gmt_files != null) {
                gmt_files = Channel.from(params.enrich.gmt_files.split(",")).collect()
            } else {
                if (params.cnmf.run_gsea || params.cnmf.run_ora) {
                    throw new Exception('No enrich.gmt files specified, but cnmf.run_gsea and or cnmf.run_ora is true')
                }
            }

            // Universe channel
            if (params.enrich.universe) {
                universe = file(params.enrich.universe)
                if (!universe.exists()) {
                    throw new Exception("Supplied id universe file does not exist")
                }
                // Custom conversion file
                universe = Channel.value(file(params.enrich.universe))
            } else {
                // In case there is no universe file
                universe = Channel.value(file("NO_UNIVERSE"))
            }
    
            if (params.cnmf.k_ignore != null) {
                gsea_in = cnmf_out.spectra_score.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                gsea_in = cnmf_out.spectra_score.map{i -> tuple("k_"+i[0], i[1])}
            }

            //----------------------------------------------------------------------------------
            // Run GSEA
            if (params.cnmf.run_gsea) {
                gsea_out = gsea("cnmf/consensus/${params.rn_runname}",
                    gsea_in,
                    gmt_files,
                    universe,
                    true)
            } else {
                gsea_out = Channel.empty()
            }
            
            //----------------------------------------------------------------------------------
            // Run ORA
            if (params.cnmf.run_ora) {
                ora_out = ora("cnmf/consensus/${params.rn_runname}",
                    gsea_in,
                    gmt_files,
                    universe,
                    true,
                    false,
                    0,
                    params.enrich.use_top)   
            } else {
                ora_out = Channel.empty()
            }

            //----------------------------------------------------------------------------------
            // Run decoupler (Progeny + collecTRI)
            if (params.cnmf.run_decoupler) {
                
                if (params.cnmf.k_ignore != null) {
                    decoupler_in = cnmf_out.spectra_score.filter { tuple ->
                        def (k, path) = tuple
                        !(k in params.cnmf.k_ignore.split(','))
                    }.map{i -> tuple("k_"+i[0], i[1])}
                } else {
                    decoupler_in = cnmf_out.spectra_score.map{i -> tuple("k_"+i[0], i[1])}
                }
                is_ensembl = (params.convert.is_ensembl_id && !params.convert.convert_gene_names) || (!params.convert.is_ensembl_id && params.convert.convert_gene_names)          
                if (is_ensembl) {
                    // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                    decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, id_linker_inv) 
                } else {
                    // This assumes you have converted to gene symbols
                    decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, file("NO_MAPPING"))
                }
            } else {
                decoupler_out = Channel.empty()
            }

            
            //----------------------------------------------------------------------------------
            // Run magma
            if (params.cnmf.run_magma) {
                
                magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)
                
                // Magma association with cnmf
                magma_cnmf_in = cnmf_out.spectra_score.filter { tuple ->
                        def (k, path) = tuple
                        !(k in params.cnmf.k_ignore.split(','))
                    }.map{i -> tuple("k_"+i[0], i[1])}
                    
                magma_assoc_in = magma_base.out.raw.combine(magma_cnmf_in)
                
                // Run regression based magma
                magma_assoc_out = magma_assoc(magma_assoc_in, true, universe, file("NO_MAPPING"))
                    
                // Concat the results in a single table
                concat_in = magma_assoc_out.out.collect().map{ list -> ["${params.rn_runname}", list]}
                magma_out = magma_concat_main("", concat_in)
                
                // Collect the results per k value as well
                magma_out_per_k = magma_concat_per_k("cnmf/consensus/${params.rn_runname}", magma_assoc_out.per_database.groupTuple())

            } else {
                magma_out = Channel.empty()
                magma_out_per_k = Channel.empty()
            }

            //----------------------------------------------------------------------------------
            // Collect enrichment results for a summary
            // Merge all the output files and calculate FDR, optionally annotate

            // For cnmf workflow this doesn't make much sense to have, so default to not annotating
            // If annotation is needed, run cnmf, then manually run_enrich later
            annot_file = Channel.value(file("NO_ANNOT"))
            
            merge_in = Channel.empty()
            .mix(gsea_out.std.map{row -> row[1]}, ora_out.std.map{row -> row[1]}, decoupler_out.std.map{row -> row[1]}, magma_out.std.map{row -> row[1]})
            .collect()
            .map{ list -> ["${params.rn_runname}", list]}
            
            concat_enrichment_results_main("cnmf/consensus/", merge_in, annot_file)

            // Per k value
            merge_in_per_k = Channel.empty()
            .mix(gsea_out.std, ora_out.std, decoupler_out.std, magma_out_per_k.std)
            .groupTuple(by:0)
            
            concat_enrichment_results_per_k("cnmf/consensus/${params.rn_runname}", merge_in_per_k, annot_file)

        } 

}