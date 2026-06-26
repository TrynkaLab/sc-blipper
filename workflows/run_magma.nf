#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; concat_enrichment_results } from "../processes/enrichment.nf"
include { preprocess_matrix } from "../processes/utils.nf"
include { prepare_magma_matrix; magma_assoc; magma_assoc as magma_assoc_permuted; magma_enrich; magma_concat;
          generate_permuted_matrix; merge_permuted_gsa; magma_permutation_stats; plot_magma_permutation;
          concat_permutation_stats } from "../processes/magma.nf"

// Subworkflows
include { perpare_enrichment } from "../subworkflows/prepare_enrichment.nf"
include { magma_base } from "../subworkflows/magma.nf"


// Main workflow
workflow magma {
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
        // Run MAGMA

        // Handle gmt and other separately
        input_matrix_split = input_matrix.branch { it ->
            gmt: it[1].extension == "gmt"
            mat: ['tsv', 'csv', 'txt'].contains(it[1].extension)
            other: false
        }

        // Prepare magma, run gene set associations
        magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)

        // Prepare input matrix once (ID mapping, universe filtering) then split for assoc and permutation
        mat_ch = input_matrix_split.mat.map { row -> tuple(row[0], file(row[1])) }
        prepared_matrix = prepare_magma_matrix(mat_ch, false, universe, file("NO_MAPPING"))
        prepared_split = prepared_matrix.multiMap { db, f ->
            for_assoc: tuple(db, f)
            for_perm:  tuple(db, f)
        }

        // Magma association with input matrix
        magma_assoc_in = magma_base.out.raw.combine(prepared_split.for_assoc)
        magma_assoc_raw = magma_assoc(magma_assoc_in)

        // Split assoc output so it can reach both magma_concat and (optionally) permutation stats
        assoc_out_split = magma_assoc_raw.out.multiMap { f ->
            concat:     f
            perm_stats: f
        }

        // Magma association with input gmt file
        magma_enrich_in = magma_base.out.raw.combine(input_matrix_split.gmt.map { row -> tuple(row[0], file(row[1])) })
        magma_enrich_out = magma_enrich(magma_enrich_in, universe, file("NO_MAPPING"))

        // Concat the results in a single table
        concat_in = assoc_out_split.concat.mix(magma_enrich_out.out).collect().map { list -> [params.rn_runname, list] }
        magma_out = magma_concat("enrich", concat_in)

        //----------------------------------------------------------------------------------
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
