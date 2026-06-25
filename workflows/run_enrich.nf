#!/usr/bin/env nextflow

// Processes
include { gsea; ora; decoupler; concat_enrichment_results; fetch_omnipath } from "../processes/enrichment.nf"
include { magma_assoc; magma_assoc as magma_assoc_permuted; magma_concat;
          generate_permuted_matrix; magma_permutation_stats; plot_magma_permutation;
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

            // Split input_matrix so it can be consumed by both real assoc and permutation
            mat_split = input_matrix
                .map { row -> tuple(row[0], file(row[1])) }
                .multiMap { db, f ->
                    real: tuple(db, f)
                    perm: tuple(db, f)
                }

            // Magma association with input matrix
            magma_assoc_in = magma_base.out.raw.combine(mat_split.real)
            magma_assoc_raw = magma_assoc(magma_assoc_in, false, universe, file("NO_MAPPING"))

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

                // Generate single permuted matrix (all N permutations per condition as separate columns)
                perm_matrix = generate_permuted_matrix(mat_split.perm, params.magma.n_permutations, params.magma.permutation_seed)

                // Run magma assoc on permuted matrix; suffix database name to keep outputs distinct
                perm_assoc_in = magma_base.out.raw.combine(
                    perm_matrix.map { db, f -> tuple("${db}_permuted", f) }
                )
                perm_assoc_out = magma_assoc_permuted(perm_assoc_in, false, universe, file("NO_MAPPING"))

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
