#!/usr/bin/env nextflow

// Subworkflows
include { fetch_id_linker } from "../subworkflows/id_linking.nf"
include { convert_and_merge } from "../subworkflows/convert_merge.nf"

// Main workflow
workflow convert {
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
        convert_and_merge(params.rn_manifest, params.rn_runname, id_linker)
                
}