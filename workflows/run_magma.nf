#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; concat_enrichment_results } from "../processes/enrichment.nf"
include { preprocess_matrix } from "../processes/utils.nf"
include { magma_assoc; magma_enrich; magma_concat;
          generate_permuted_matrix; magma_permutation_stats; plot_magma_permutation } from "../processes/magma.nf"

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

        // Split mat channel so it can be consumed by both real assoc and permutation
        mat_split = input_matrix_split.mat
            .map { row -> tuple(row[0], file(row[1])) }
            .multiMap { db, f ->
                real: tuple(db, f)
                perm: tuple(db, f)
            }

        // Prepare magma, run gene set associations
        magma_base(params.magma, params.convert, params.rn_ensembl_version, ensembl_reference)

        // Magma association with input matrix
        magma_assoc_in = magma_base.out.raw.combine(mat_split.real)
        magma_assoc_raw = magma_assoc(magma_assoc_in, false, universe, file("NO_MAPPING"))

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

            // Generate single permuted matrix (all N permutations per condition as separate columns)
            perm_matrix = generate_permuted_matrix(mat_split.perm, params.magma.n_permutations)

            // Run magma assoc on permuted matrix; suffix database name to keep outputs distinct
            perm_assoc_in = magma_base.out.raw.combine(
                perm_matrix.map { db, f -> tuple("${db}_permuted", f) }
            )
            perm_assoc_out = magma_assoc(perm_assoc_in, false, universe, file("NO_MAPPING"))

            // Parse (trait, database) from gsa.out filenames for channel joining
            real_gsa = assoc_out_split.perm_stats.map { f ->
                def s = f.name.replace('.gsa.out', '')
                def i = s.indexOf('__')
                tuple(s.substring(0, i), s.substring(i + 2), f)
            }

            perm_gsa = perm_assoc_out.out.map { f ->
                def s = f.name.replace('.gsa.out', '')
                def i = s.indexOf('__')
                def trait = s.substring(0, i)
                def db_perm = s.substring(i + 2)
                def db = db_perm.substring(0, db_perm.lastIndexOf('_permuted'))
                tuple(trait, db, f)
            }

            // Split perm_gsa for use in both stats and plot processes
            perm_gsa_split = perm_gsa.multiMap { trait, db, f ->
                stats: tuple(trait, db, f)
                plot:  tuple(trait, db, f)
            }

            // Compute permutation statistics
            stats_in = real_gsa.join(perm_gsa_split.stats, by: [0, 1])
            perm_stats_out = magma_permutation_stats(stats_in)

            // Plot null distributions
            plot_in = perm_stats_out.join(perm_gsa_split.plot, by: [0, 1])
            plot_magma_permutation(plot_in)
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
