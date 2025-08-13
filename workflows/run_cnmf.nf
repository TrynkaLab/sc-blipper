#!/usr/bin/env nextflow

include { cnmf_pre_process; cnmf_prepare; cnmf_factorize; cnmf_combine; cnmf_kselection; cnmf_consensus } from "../processes/cnmf.nf"
include { gsea; ora; decoupler } from "../processes/enrichment.nf"
include { magma_assoc } from "../processes/magma.nf"

// Sub workflows
include { fetch_id_linker } from "../workflows/id_linking.nf"
include { convert_and_merge } from "../workflows/convert_merge.nf"
include { magma_base } from "../workflows/magma.nf"

// Main workflow
workflow cnmf {
    //TODO: make gene linker file ids unique
    main:

        //------------------------------------------------------------
        // Id linking and ensembl reference
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        
        // Add the output in the global space for convenience
        id_linker = fetch_id_linker.out.id_linker
        id_linker_inv = fetch_id_linker.out.id_linker_inv
        ensembl_reference = fetch_id_linker.out.ensembl_reference

        //------------------------------------------------------------        
        convert_and_merge(params.rn_manifest, params.rn_runname)
        
        merge_out = convert_and_merge.out.merge_out
        
        //--------------------------------------------------------
        // Optionally run cNMF preprocessing, which creates harmony corrected counts
        if (params.cnmf.preprocess) {
            // Preprocess the cnmf file
            cnmf_preprocess = cnmf_pre_process(merge_out)
            
            // Channel with run name and .h5ad file with counts
            cnmf_in = cnmf_preprocess.preproccessed
        } else {
            // Otherwise use the merged file, in case there is one batch,
            // just uses the other file
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
         .map { v -> (params.cnmf.k.split(",")).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out)
        // This is the end of the cnmf processing

        //--------------------------------------------------------
        // Optional post-proccessing of cnmf
        //--------------------------------------------------------

        //----------------------------------------------------------------------------------
        // Run gene set enrichment
        if (params.cnmf.run_enrichment && params.enrich.gmt_files != null) {
            
            // Channel with gmt files
            gmt_files = Channel.from(params.enrich.gmt_files.split(",")).collect()

            // Universe channel
            if (params.enrich.universe) {
                universe = file(params.enrich.universe)
                if (!id_linker.exists()) {
                    throw new Exception("Supplied id universe file does not exist")
                }
                // Custom conversion file
                universe = Channel.value(file(params.enrich.universe))
            } else {
                // In case there is no universe file
                universe = Channel.value(file("NO_UNIVERSE"))
            }
    
            if (params.cnmf.k_ignore != null) {
                gsea_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                gsea_in = cnmf_out.spectra_k.map{i -> tuple("k_"+i[0], i[1])}
            }

            // Run GSEA
            gsea_out = gsea("cnmf/consensus/${params.rn_runname}",
                            gsea_in,
                            gmt_files,
                            universe,
                            true)
                            
            // Run ORA
            ora_out = ora("cnmf/consensus/${params.rn_runname}",
                gsea_in,
                gmt_files,
                universe,
                true,
                false,
                0,
                params.enrich.use_top)
        }
        
        //----------------------------------------------------------------------------------
        // Run decoupler (Progeny + collecTRI)
        if (params.cnmf.run_decoupler) {
            
            if (params.cnmf.k_ignore != null) {
                decoupler_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                decoupler_in = cnmf_out.spectra_k.map{i -> tuple("k_"+i[0], i[1])}
            }
            
            if (params.convert.invert_linker) {
                // This assumes you have converted to gene symbols
                decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, file("NO_MAPPING"))
            } else {
                // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, id_linker_inv) 
            }
        }
        
        //----------------------------------------------------------------------------------
        // Run magma
        if (params.cnmf.run_magma) {
            
            magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)
            
            // Magma association with cnmf
            magma_cnmf_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
                
            magma_assoc_in = magma_base.out.raw.combine(magma_cnmf_in)
            
            // Run regression based magma
            magma_assoc_out = magma_assoc(magma_assoc_in, true)
         
            // Concat the results in a single table
            magma_concat(params.rn_runname, magma_assoc_out.out.collect())
            
        }
        
        
        // Collect enrichment results for a summary

}