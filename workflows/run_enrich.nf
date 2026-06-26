#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; concat_enrichment_results; fetch_omnipath } from "../processes/enrichment.nf"
include { prepare_magma_matrix; magma_assoc; magma_assoc as magma_assoc_permuted; magma_concat;
          generate_permuted_matrix; merge_permuted_gsa; magma_permutation_stats; plot_magma_permutation;
          concat_permutation_stats } from "../processes/magma.nf"

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
            def use_top = params.enrich.use_top != null ? params.enrich.use_top : []

            ora_out = ora("enrich/",
                input_matrix,
                gmt_files,
                universe,
                false,
                params.enrich.absolute,
                params.enrich.threshold,
                params.enrich.threshold_invert,
                use_top
            )
        } else {
            ora_out = [:]
            ora_out.std = Channel.empty()
        }

        //----------------------------------------------------------------------------------
        // Run decoupler (Progeny + collecTRI)
        if (params.enrich.run_decoupler) {

            // Fetch the reference files
            if (params.enrich.progeny_file != null && params.enrich.collectri_file != null) {
                progeny_network = Channel.value(file(params.enrich.progeny_file))
                collectri_network = Channel.value(file(params.enrich.collectri_file))
            } else if (params.enrich.progeny_file != null && params.enrich.collectri_file == null || params.enrich.progeny_file == null && params.enrich.collectri_file != null) {
                error "Must specify both progeny and collectri files, or neither to fetch from omnipath"
            } else {
                omnipath = fetch_omnipath()
                progeny_network = omnipath.progeny
                collectri_network = omnipath.collectri
            }

            is_ensembl = params.convert.output_namespace == "ensembl"

            if (is_ensembl) {
                // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                decoupler_out = decoupler("enrich/", input_matrix, false, id_linker_inv, progeny_network, collectri_network)
            } else {
                // This assumes you have converted to gene symbols
                decoupler_out = decoupler("enrich/", input_matrix, false, file("NO_MAPPING"), progeny_network, collectri_network)
            }
        } else {
            decoupler_out = [:]
            decoupler_out.std = Channel.empty()
        }

        //----------------------------------------------------------------------------------
        // Run MAGMA
        if (params.enrich.run_magma) {
            magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)

            // Prepare input matrix once (ID mapping, universe filtering) then split for assoc and permutation
            mat_ch = input_matrix.map { row -> tuple(row[0], file(row[1])) }
            prepared_matrix = prepare_magma_matrix(mat_ch, false, universe, file("NO_MAPPING"))
            prepared_split = prepared_matrix.multiMap { db, f ->
                for_assoc: tuple(db, f)
                for_perm:  tuple(db, f)
            }

            // Magma association with input matrix
            magma_assoc_in = magma_base.out.raw.combine(prepared_split.for_assoc)
            magma_assoc_raw = magma_assoc(magma_assoc_in)

            // Split assoc output for concat and (optionally) permutation stats
            assoc_out_split = magma_assoc_raw.out.multiMap { f ->
                concat:     f
                perm_stats: f
            }

            // Concat the results in a single table
            concat_in = assoc_out_split.concat.collect().map { list -> [params.rn_runname, list] }
            magma_out = magma_concat("enrich", concat_in)

            // Permutation testing (set params.magma.run_permutation_test = true to enable)
            if (params.magma.run_permutation_test) {

                // Compute block ranges from permutation_block_ratio
                def n         = params.magma.n_permutations
                def blockSize = Math.max(1, (int) Math.ceil(params.magma.permutation_block_ratio * n))
                def blocks    = []
                def blkStart  = 1
                def blkIdx    = 0
                while (blkStart <= n) {
                    def blkEnd = Math.min(blkStart + blockSize - 1, n)
                    blocks << [blkIdx, blkStart, blkEnd]
                    blkStart = blkEnd + 1
                    blkIdx++
                }
                def perm_blocks_ch = Channel.fromList(blocks)

                // Cross each (database, matrix) with each block spec, then generate per-block permuted matrices
                perm_matrix_in = prepared_split.for_perm.combine(perm_blocks_ch)
                perm_matrix = generate_permuted_matrix(perm_matrix_in, params.magma.n_permutations, params.magma.permutation_seed)

                // Run magma assoc per block; tag db name with block index to keep outputs distinct
                perm_assoc_in = magma_base.out.raw.combine(
                    perm_matrix.matrix.map { db, bidx, f -> tuple("${db}_permuted_block${bidx}", f) }
                )
                perm_assoc_out = magma_assoc_permuted(perm_assoc_in)

                // Parse (trait, db) from filename — strip _permuted_block{N} suffix
                perm_gsa_blocks = perm_assoc_out.out.map { f ->
                    def s      = f.name.replace('.gsa.out', '')
                    def i      = s.indexOf('__')
                    def trait  = s.substring(0, i)
                    def db_raw = s.substring(i + 2)
                    def db     = db_raw.replaceAll(/_permuted_block\d+$/, '')
                    tuple(trait, db, f)
                }

                // Group all block outputs for the same (trait, db), then merge into a single gsa.out
                merged_perm_gsa = merge_permuted_gsa(perm_gsa_blocks.groupTuple(by: [0, 1]))

                // Parse (trait, database) from real gsa.out filenames for channel joining
                real_gsa = assoc_out_split.perm_stats.map { f ->
                    def s = f.name.replace('.gsa.out', '')
                    def i = s.indexOf('__')
                    tuple(s.substring(0, i), s.substring(i + 2), f)
                }

                // Split merged perm_gsa for use in both stats and plot processes
                perm_gsa_split = merged_perm_gsa.multiMap { trait, db, f ->
                    stats: tuple(trait, db, f)
                    plot:  tuple(trait, db, f)
                }

                // Compute permutation statistics
                stats_in = real_gsa.join(perm_gsa_split.stats, by: [0, 1])
                perm_stats_out = magma_permutation_stats(stats_in)

                // Split stats output for plot and concat
                perm_stats_split = perm_stats_out.multiMap { trait, db, f ->
                    for_plot:   tuple(trait, db, f)
                    for_concat: f
                }

                // Plot null distributions
                plot_in = perm_stats_split.for_plot.join(perm_gsa_split.plot, by: [0, 1])
                plot_magma_permutation(plot_in)

                // Concat all permutation stats into a single file
                concat_perm_in = perm_stats_split.for_concat.collect().map { files -> tuple(params.rn_runname, files) }
                concat_permutation_stats(concat_perm_in)
            }

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
