#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; preprocess_matrix; concat_enrichment_results } from "../processes/enrichment.nf"
include { magma_assoc; magma_concat } from "../processes/magma.nf"

// Subworkflows
include { perpare_enrichment } from "../subworkflows/prepare_enrichment.nf"
include { magma_base } from "../subworkflows/magma.nf"


// Main workflow
workflow enrich {
    main:
    
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
        // Run GSEA
        if (params.enrich.run_gsea) {
            gsea_out = gsea("enrich/",
                input_matrix,
                gmt_files,
                universe,
                false)             
        } else {
            gsea_out = [:]
            gsea_out.std = Channel.empty()
        }
        
        //----------------------------------------------------------------------------------
        // Run ORA
        if (params.enrich.run_ora) {
            ora_out = ora("enrich/",
                input_matrix,
                gmt_files,
                universe,
                false,
                params.enrich.absolute,
                params.enrich.threshold,
                params.enrich.use_top)
        } else {
            ora_out = [:]
            ora_out.std = Channel.empty()
        }
        
        //----------------------------------------------------------------------------------
        // Run decoupler (Progeny + collecTRI)
        if (params.enrich.run_decoupler) {  
            is_ensembl = (params.convert.is_ensembl_id && !params.convert.convert_gene_names) || (!params.convert.is_ensembl_id && params.convert.convert_gene_names)          
            if (is_ensembl) {
                // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                decoupler_out = decoupler("enrich/", input_matrix, false, id_linker_inv) 
            } else {
                // This assumes you have converted to gene symbols
                decoupler_out = decoupler("enrich/", input_matrix, false, file("NO_MAPPING"))
            }
        } else {
            decoupler_out = [:]
            decoupler_out.std = Channel.empty()
        }
        
        //----------------------------------------------------------------------------------
        // Run MAGMA
        if (params.enrich.run_magma) {
            magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)
            
            // Magma association with input matrix
            magma_assoc_in = magma_base.out.raw.combine(input_matrix.map{row -> tuple(row[0], file(row[1]))})
                        
            // Run regression based magma
            magma_assoc_out = magma_assoc(magma_assoc_in, false, universe, file("NO_MAPPING"))
         
            // Concat the results in a single table
            concat_in = magma_assoc_out.out.collect().map{ list -> [params.rn_runname, list]}
            magma_out = magma_concat("enrich", concat_in)
        } else {
            magma_out = [:]
            magma_out.std = Channel.empty()
        }
        
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
        .mix(gsea_out.std,
             ora_out.std,
             decoupler_out.std.flatMap { id, files -> files.collect { file -> tuple(id, file) }},
             magma_out.std)
        .groupTuple(by:0)
    
        // merge_in.view()
        concat_enrichment_results("enrich", merge_in, annot_file)
        
        
        
        
        
    
}
