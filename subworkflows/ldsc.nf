#!/usr/bin/env nextflow
include { munge_sumstats; ensembl_to_ldsc_geneloc; make_annotation_and_ldscores; make_annotation_and_ldscores_per_chr; collect_ldscores; ldsc_stratified; aggregate_ldsc_results } from '../processes/ldsc.nf'
include { binarize_matrix; extract_conditions; preprocess_matrix } from '../processes/utils.nf'
include { fetch_id_linker } from '../subworkflows/id_linking.nf'


// Partitioned (stratified) LD score regression subworkflow.
//
// Produces enrichment statistics showing what fraction of SNP heritability
// (Prop._h2 in the output) is captured by each gene set column.
//
// MANIFEST FORMAT (tab-separated, one row per phenotype):
//   name              – unique phenotype identifier
//   path              – path to the raw summary statistics file
//   variant_col       – column name for the SNP / variant ID
//   effect_allele_col – column name for the effect allele (A1)
//   other_allele_col  – column name for the other allele (A2)
//   sign_col          – column name for the signed statistic passed to --signed-sumstats
//                       e.g. beta column, OR column, or Z score column
//   sign_null         – null value for the signed statistic (0 for beta/Z, 1 for OR)
//   p_col             – column name for the p-value
//   n                 – fixed total sample size;  ignored if n_case and n_control are both set
//   n_case            – fixed N cases;            use "NA" if not a case/control study
//   n_control         – fixed N controls;         use "NA" if not a case/control study
//
// GENE MATRIX FORMAT (tab-separated):
//   First column: gene identifiers (row names).
//   Remaining columns: one column per condition.  Two modes:
//     Binary matrix  (params.ldsc.binarize_top = null):
//       values must be 0/1; 1 = gene in set, 0 = background.
//     Numeric matrix (params.ldsc.binarize_top set):
//       any numeric scores; top N genes per column are selected as the gene
//       set.  Set params.ldsc.binarize_absolute to take |scores| first, and
//       params.ldsc.binarize_ascending to select the lowest scores instead.
//
// REFERENCE FOLDER LAYOUT:
//   hapmap3/w_hm3.snplist     – HapMap3 SNP list for munge_sumstats
//   hapmap3/hm3_no_MHC.list   – SNP IDs to restrict LD score computation (MHC excluded)
//   plink/{plink_prefix}{1..22}.{bed,bim,fam}
//   baselineLD_v2.2/baselineLD.{1..22}.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}
//   weights/weights.hm3_noMHC.{1..22}.l2.ldscore.gz
//   frq/{frq_prefix}{1..22}.frq
//
// CONDA / CONTAINER:
//   Requires the CBIIT Python 3.9 LDSC environment (github.com/CBIIT/ldsc).
//   Set params.ldsc.conda to the path of that environment.
//
// PARAMS (set in nextflow.config under params.ldsc):
//   label          – process label for resource allocation  (e.g. "normal")
//   label_high     – label for the annotation + LD score step (e.g. "normal_plus")
//   conda          – path to the CBIIT ldsc conda environment
//   window_size    – gene window in kb for make_annot.py     (default 100)
//   ld_wind_cm     – LD window in cM for --l2 step           (default 1)
//   plink_prefix   – filename prefix for plink files in plink/  (default "chr")
//   frq_prefix     – filename prefix for frq files in frq/      (default "frq.")
//
//   Override reference paths (optional; default derived from reference_dir):
//   plink_dir      – full prefix up to chromosome number for plink files
//                    (default: reference_dir/plink/{plink_prefix})
//   hm3_no_mhc     – path to hm3_no_MHC.list for LD score computation
//                    (default: reference_dir/hapmap3/hm3_no_MHC.list)
//   w_hm3_snplist  – path to w_hm3.snplist for munge_sumstats
//                    (default: reference_dir/hapmap3/w_hm3.snplist)
//   baseline_ld_chr – chr-split prefix for baselineLD ref-ld-chr
//                    (default: reference_dir/baselineLD_v2.2/baselineLD.)
//   weights_chr    – chr-split prefix for weights w-ld-chr
//                    (default: reference_dir/weights/weights.hm3_noMHC.)
//   frq_chr        – chr-split prefix for frqfile-chr
//                    (default: reference_dir/frq/{frq_prefix})

workflow ldsc {

    take:
        manifest       // string path: summary statistics manifest TSV
        gene_matrix    // string path: genes x conditions matrix (0/1 coded)
        reference_dir  // string path: reference folder (see layout above)

    main:

        //------------------------------------------------------------
        // Validate required parameters
        //------------------------------------------------------------
        if (params.ldsc == null) {
            throw new Exception("params.ldsc block is not configured. Add it to nextflow.config.")
        }
        if (params.ldsc.conda == null) {
            throw new Exception("params.ldsc.conda must point to the CBIIT ldsc conda environment.")
        }
        if (params.ldsc.window_size == null) {
            throw new Exception("params.ldsc.window_size must be set (gene window in kb, e.g. 100).")
        }

        //------------------------------------------------------------
        // Munge summary statistics
        //------------------------------------------------------------
        manifest_ch = Channel.fromPath(manifest)
            .splitCsv(header: true, sep: "\t")
            .map { row -> tuple(
                row.name,
                file(row.path),
                row.variant_col        ?: "SNP",
                row.effect_allele_col  ?: "A1",
                row.other_allele_col   ?: "A2",
                row.sign_col,
                row.sign_null,
                row.p_col,
                row.n                  ?: "NA",
                row.n_case             ?: "NA",
                row.n_control          ?: "NA"
            )}

        munged_ch = munge_sumstats(manifest_ch, reference_dir).munged

        //------------------------------------------------------------
        // Fetch Ensembl reference — supplies both gene coords and id_linker
        //------------------------------------------------------------
        is_ensembl = params.convert.output_namespace == "ensembl"
        fetch_id_linker(params.rn_ensembl_version, params.convert)
        gene_coords = ensembl_to_ldsc_geneloc(fetch_id_linker.out.ensembl_reference, is_ensembl).gene_coords


        //------------------------------------------------------------
        // Preprocess matrix: standardise gene IDs via id_linker
        //------------------------------------------------------------
        // TODO: Double check the gene id conversion works propely
        prep_in = Channel.value(tuple(params.rn_runname, gene_matrix, params.enrich.transpose))
            .combine(fetch_id_linker.out.id_linker)
            .combine(Channel.value(file("NO_UNIVERSE")))
            .map { id, path, transpose, id_lnk, universe ->
                tuple(id, file(path), transpose, id_lnk, universe)
            }
            
        preprocessed_ch = preprocess_matrix("ldsc/", prep_in, file("NO_SUBSET")).matrix
            .map { id, matrix -> matrix }

        //------------------------------------------------------------
        // Resolve window_size and binarize_top as lists
        //------------------------------------------------------------
        def ws_list = params.ldsc.window_size instanceof List
            ? params.ldsc.window_size
            : [params.ldsc.window_size]
        def bt_list = params.ldsc.binarize_top != null
            ? (params.ldsc.binarize_top instanceof List ? params.ldsc.binarize_top : [params.ldsc.binarize_top])
            : [null]

        // Use .first() so the single preprocessed matrix can be combined many times
        preprocessed_val = preprocessed_ch.first()

        //------------------------------------------------------------
        // Binarize for each top value (or pass through with null tag)
        //------------------------------------------------------------
        if (params.ldsc.binarize_top != null) {
            binarize_input_ch = Channel.from(bt_list).combine(preprocessed_val)
            binarize_input_ch.multiMap { bt, m -> tops: bt; matrices: m }.set { bsplit }

            binarized_ch = binarize_matrix(
                "ldsc/",
                params.rn_runname,
                bsplit.matrices,
                bsplit.tops,
                params.ldsc.binarize_absolute,
                params.ldsc.binarize_ascending,
            ).matrix
                .map { m ->
                    // Extract top value as string to avoid null-check issues downstream
                    def top = (m.name =~ /__top(\d+)\.tsv$/)[0][1]
                    tuple(top, m)
                }
        } else {
            binarized_ch = preprocessed_val.map { m -> tuple("NA", m) }
        }

        //------------------------------------------------------------
        // Split each binarized matrix into per-condition gene set files
        //------------------------------------------------------------
        // extract_conditions carries the binarize_top tag through
        gene_sets_ch = extract_conditions(binarized_ch).gene_sets
            .flatMap { bt, files ->
                (files instanceof List ? files : [files])
                    .collect { f -> tuple(bt, f.name.replace(".geneset.txt", ""), f) }
            }
        // gene_sets_ch: (binarize_top, condition_name, gene_set_file)

        // Cross with every window_size
        annotate_input_ch = gene_sets_ch
            .combine(Channel.from(ws_list))
            .map { bt, cond, gs, ws -> tuple(ws, bt, cond, gs) }
        // annotate_input_ch: (window_size, binarize_top, condition_name, gene_set_file)

        //------------------------------------------------------------
        // Create annotation files and compute LD scores
        //------------------------------------------------------------
        if (params.ldsc.calculate_ldscores_per_chr) {
            per_chr_input_ch = annotate_input_ch
                .combine(Channel.from(1..22))
                .map { ws, bt, cond, gs, chr -> tuple(cond, gs, chr, ws, bt) }

            ldscores_ch = make_annotation_and_ldscores_per_chr(
                    per_chr_input_ch, reference_dir, gene_coords
                ).ldscores
                .groupTuple(by: [0, 1, 2])
                .map { cond, ws, bt, chrs, files -> tuple(cond, ws, bt, files.flatten()) }
            ldscores_ch = collect_ldscores(ldscores_ch).ldscores
        } else {
            annot_input_ch = annotate_input_ch.map { ws, bt, cond, gs -> tuple(cond, gs, ws, bt) }
            ldscores_ch = make_annotation_and_ldscores(annot_input_ch, reference_dir, gene_coords).ldscores
        }
        // ldscores_ch: (condition_name, window_size, binarize_top, ldscores_dir)

        //------------------------------------------------------------
        // Run partitioned LDSC for every phenotype x condition x ws x bt
        //------------------------------------------------------------
        ldsc_input_ch = munged_ch
            .combine(ldscores_ch)
            .map { pheno_id, munged, condition_name, ws, bt, ldscores ->
                tuple(pheno_id, munged, condition_name, ldscores, ws, bt)
            }

        results_ch = ldsc_stratified(ldsc_input_ch, reference_dir).results
        // results_ch: (pheno_id, condition_name, window_size, binarize_top, results_file)

        //------------------------------------------------------------
        // Collect all .results files and aggregate into one table
        //------------------------------------------------------------
        results_files_ch = results_ch
            .map { pheno_id, condition_name, ws, bt, results_file -> results_file }
            .collect()

        aggregated = aggregate_ldsc_results(results_files_ch).aggregated

    emit:
        aggregated = aggregated   // single TSV with all phenotypes, conditions, window_sizes, binarize_tops
        results    = results_ch   // tuple(pheno_id, condition_name, window_size, binarize_top, results_file)

}
