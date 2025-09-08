#!/usr/bin/env nextflow
include { fetch_gene_id_reference; invert_id_link } from "../processes/utils.nf"
include { fetch_id_linker } from "../subworkflows/id_linking.nf"

workflow perpare_enrichment {
    take:
        params
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
        if (params.enrich.gmt_files != null) {
            gmt_files = Channel.from(params.enrich.gmt_files.split(",")).collect()
        } else {
            gmt_files = Channel.empty()
        }

        // Universe channel
        if (params.enrich.universe) {
            universe = file(params.enrich.universe)
            if (!universe.exists()) {
                throw new Exception("Supplied id universe file does not exist")
            }
            
            if (params.convert.convert_gene_names) {
                log.warn("Universe is not converted to output name space, so make sure it is matched to what you are converting to")
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
        //input_matrix = Channel.value(tuple(params.rn_runname, params.enrich.input_matrix))
        if (params.convert.convert_gene_names) {
            converter = id_linker
        } else {
            converter = Channel.value(file("NO_MAPPING"))
        }
        
        
        if (params.enrich.input_matrix == null) {
            throw new Exception("params.enrich.input_matrix cannot be null")
        } else {
            
            if (params.enrich.input_matrix.endsWith(".gmt")) {
                input_matrix = Channel.value(tuple(params.rn_runname, file(params.enrich.input_matrix)))
                
            } else {
                prep_in = Channel.value(tuple(params.rn_runname, params.enrich.input_matrix, params.enrich.transpose))
                .combine(converter)
                .combine(universe)
                .map{row -> tuple(row[0], file(row[1]), row[2], file(row[3]), file(row[4]))}

                input_matrix = preprocess_matrix("enrich/", prep_in).matrix
            }
        }
        

        
    emit:
        input_matrix = input_matrix
        converter = converter
        universe = universe
        gmt_files = gmt_files
        id_linker = id_linker
        id_linker_inv = id_linker_inv
        ensembl_reference = ensembl_reference
}