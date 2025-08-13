#!/usr/bin/env nextflow
include { fetch_gene_id_reference; invert_id_link } from "../processes/utils.nf"

workflow fetch_id_linker {

    take:
        version
        convert_params

    main:
        // Always fetch this even if not used
        ensembl_reference = fetch_gene_id_reference(version)

        if (convert_params.convert_gene_names) {
            if (convert_params.id_linker) {
                id_file = file(convert_params.id_linker)
                if (!id_file.exists()) {
                    error "Supplied id linker file does not exist: ${convert_params.id_linker}"
                }
                id_linker     = Channel.value(id_file)
                id_linker_inv = invert_id_link(id_linker)

            } else {
                if (convert_params.invert_linker) {
                    // ENSEMBL → gene name
                    id_linker     = ensembl_reference.ensembl_to_name
                    id_linker_inv = ensembl_reference.name_to_ensembl
                } else {
                    // gene name → ENSEMBL
                    id_linker     = ensembl_reference.name_to_ensembl
                    id_linker_inv = ensembl_reference.ensembl_to_name
                }
            }
        }
        else {
            // No mapping
            dummy_file = file("NO_MAPPING")
            id_linker     = Channel.value(dummy_file)
            id_linker_inv = Channel.value(dummy_file)
        }

    emit:
        ensembl_reference = ensembl_reference.ensembl
        id_linker = id_linker
        id_linker_inv = id_linker_inv
}
