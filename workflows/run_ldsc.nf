#!/usr/bin/env nextflow

// Subworkflows
include { ldsc as ldsc_base } from "../subworkflows/ldsc.nf"


// Main workflow
workflow ldsc {
    main:

        //----------------------------------------------------------------------------------
        // Validate inputs
        if (params.ldsc.manifest_sumstats == null) {
            throw new Exception("params.ldsc.manifest_sumstats must be supplied")
        }
        if (params.ldsc.gene_matrix == null) {
            throw new Exception("params.ldsc.gene_matrix must be supplied")
        }
        if (params.ldsc.reference_dir == null) {
            throw new Exception("params.ldsc.reference_dir must be supplied")
        }

        //----------------------------------------------------------------------------------
        // Run LDSC
        ldsc_base(
            params.ldsc.manifest_sumstats,
            params.ldsc.gene_matrix,
            params.ldsc.reference_dir
        )
}
