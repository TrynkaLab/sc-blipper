#!/usr/bin/env nextflow

// Subworkflows
include { ldsc as ldsc_base } from "../subworkflows/ldsc.nf"
include { perpare_enrichment } from "../subworkflows/prepare_enrichment.nf"


// Main workflow
workflow ldsc {
    main:

        //----------------------------------------------------------------------------------
        // Validate inputs
        if (params.ldsc.manifest_sumstats == null) {
            throw new Exception("params.ldsc.manifest_sumstats must be supplied")
        }
        if (params.ldsc.reference_dir == null) {
            throw new Exception("params.ldsc.reference_dir must be supplied")
        }

        //----------------------------------------------------------------------------------
        // Prepare enrichment matrix
        perpare_enrichment(params)

        //----------------------------------------------------------------------------------
        // Run LDSC
        ldsc_base(
            params.ldsc.manifest_sumstats,
            perpare_enrichment.out.input_matrix,
            perpare_enrichment.out.ensembl_reference,
            params.ldsc.reference_dir
        )
}
