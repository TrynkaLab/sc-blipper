#!/usr/bin/env nextflow

// Subworkflows
include { fetch_id_linker } from "../subworkflows/id_linking.nf"
include { convert_and_merge } from "../subworkflows/convert_merge.nf"
include { cnmf_pre_process } from "../processes/cnmf.nf"

// Main workflow
workflow convert {
    //TODO: make gene linker file ids unique
    main:

        //------------------------------------------------------------
        // Id linking and ensembl reference
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        
        //------------------------------------------------------------        
        convert_and_merge(params.rn_manifest,
                        params.rn_runname,
                        fetch_id_linker.out.id_linker,
                        params.convert.subset_genes,
                        fetch_id_linker.out.biotype_ensembl,
                        fetch_id_linker.out.biotype_gene_name,
                        params.convert.output_namespace)
                        
        //------------------------------------------------------------
        // Optionally correct batch effects using harmony
        if (params.convert.harmony_vars != null) {
            cnmf_pre_process(convert_and_merge.out.merge_out)
        }
}