#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler } from "../processes/enrichment.nf"
include { magma_assoc } from "../processes/magma.nf"

// Subworkflows
include { fetch_id_linker } from "../subworkflows/id_linking.nf"
include { magma_base } from "../subworkflows/magma.nf"


// Main workflow
workflow enrich {
    main:
    
        if (params.enrich.input_matrix == null) {
            throw new Exception("--enrich.input_matrix must be supplied")
        }
        
        if (params.rn_runname == null) {
            throw new Exception("--rn_runname must be supplied")
        }    
        
        //------------------------------------------------------------
        // Id linking and ensembl reference
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        
        // Add the output in the global space for convenience
        id_linker = fetch_id_linker.out.id_linker
        id_linker_inv = fetch_id_linker.out.id_linker_inv
        ensembl_reference = fetch_id_linker.out.ensembl_reference
        
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
        
        // Input matrix
        // TODO make a manfifest and make these configurable
        // <name> <matrix> <transpose> <absolute> GSEA | ora_top | ora_thresh | decoupler | magma
        // foo    bar.tsv  T            T         T       F           T            T         F
        input_matrix = Channel.value(tuple(params.rn_runname, params.enrich.input_matrix))
        
        //----------------------------------------------------------------------------------
        // Run GSEA
        if (params.enrich.run_gsea) {
            gsea_out = gsea("enrich/${params.rn_runname}",
                gsea_in,
                gmt_files,
                universe,
                params.enrich.transpose)             
        }

        //----------------------------------------------------------------------------------
        // Run ORA
        if (params.enrich.run_gsea) {
            ora_out = ora("enrich/${params.rn_runname}",
                gsea_in,
                gmt_files,
                universe,
                params.enrich.transpose,
                params.enrich.absolute,
                params.enrich.threshold,
                params.enrich.use_top)
        }
        
        //----------------------------------------------------------------------------------
        // Run decoupler (Progeny + collecTRI)
        if (params.enrich.run_decoupler) {            
            if (params.convert.invert_linker) {
                // This assumes you have converted to gene symbols
                decoupler_out = decoupler("enrich/${params.rn_runname}", decoupler_in, params.enrich.transpose, file("NO_MAPPING"))
            } else {
                // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                decoupler_out = decoupler("enrich/${params.rn_runname}", decoupler_in, params.enrich.transpose, id_linker_inv) 
            }
        }
        
        //----------------------------------------------------------------------------------
        // Run MAGMA
        if (params.enrich.run_magma) {
            magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)
            
            // Magma association with input matrix
            magma_assoc_in = magma_base.out.raw.combine(input_matrix)
            
            // Run regression based magma
            magma_assoc_out = magma_assoc(magma_assoc_in, params.enrich.transpose)
         
            // Concat the results in a single table
            magma_concat(params.rn_runname, magma_assoc_out.out.collect())
        
        }
    
}
