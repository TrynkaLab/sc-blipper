#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; preprocess_matrix; concat_enrichment_results } from "../processes/enrichment.nf"
include { magma_assoc; magma_enrich; magma_concat } from "../processes/magma.nf"

// Subworkflows
include { perpare_enrichment } from "../subworkflows/prepare_enrichment.nf"
include { magma_base } from "../subworkflows/magma.nf"


// Main workflow
workflow magma {
    main:
    
        if (params.convert.convert_gene_names) {
            throw new Exception("Currently magma workflow does not support gene name conversion as converting gene names in gmt is not implemented")
        }
        
        //----------------------------------------------------------------------------------
        // Prepare the enrichment
        enrich = perpare_enrichment(params)
        input_matrix = enrich.input_matrix
        converter = enrich.converter
        universe = enrich.universe
        gmt_files = enrich.gmt_files
        id_linker = enrich.id_linker
        id_linker_inv = enrich.id_linker_inv
        ensembl_reference = enrich.ensembl_reference

        //----------------------------------------------------------------------------------
        // Run MAGMA
      
        // Handle gmt and other seperately
        input_matrix_split = input_matrix.branch{ it ->
            gmt: it[1].extension == "gmt"
            mat: ['tsv', 'csv', 'txt'].contains(it[1].extension)
            other: false
        }
        
        // Prepare magma, run gene set associations 
        magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)
            
        // Magma association with input matrix
        magma_assoc_in = magma_base.out.raw.combine(input_matrix_split.mat.map{row -> tuple(row[0], file(row[1]))})
        magma_assoc_out = magma_assoc(magma_assoc_in, false, universe, file("NO_MAPPING"))
        

        // Magma association with input gmt file
        magma_enrich_in = magma_base.out.raw.combine(input_matrix_split.gmt.map{row -> tuple(row[0], file(row[1]))})
        magma_enrich_out = magma_enrich(magma_enrich_in, universe, file("NO_MAPPING"))
        
        // Concat the results in a single table
        concat_in = magma_assoc_out.out.mix(magma_enrich_out.out).collect().map{list -> [params.rn_runname, list]}
        magma_out = magma_concat("enrich", concat_in)

        //----------------------------------------------------------------------------------
        // Merge all the output files and calculate FDR, optionally annotate
        if (params.enrich.annotate != null) {
            annot_file = Channel.value(file(params.enrich.annotate))
        } else {
            annot_file = Channel.value(file("NO_ANNOT"))
        }
        
        // Create one channel with all the files to merge. Progeny outputs multiple files, so first 
        // make each file into its own [id, file] tuple
        merge_in = Channel.empty()
        .mix(magma_out.std)
        .groupTuple(by:0)
    
        // merge_in.view()
        concat_enrichment_results("enrich", merge_in, annot_file)
}
