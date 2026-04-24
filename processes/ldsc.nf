#!/usr/bin/env nextflow


// Convert Ensembl BioMart reference to the gene coordinate file expected by
// make_annot.py: GENE  CHR  START  END (tab-separated, with header).
// Filters to protein-coding genes on autosomes + chrX and deduplicates by GENE.
process ensembl_to_ldsc_geneloc {
    label "tiny"

    container params.ldsc.container
    conda params.rn_conda

    publishDir "$params.rn_publish_dir/reference/ldsc", mode: 'symlink'

    input:
        path(ensembl_file)
        val(use_ensembl_id)

    output:
        path("gene_coords.txt"), emit: gene_coords

    script:
        def gene_col = use_ensembl_id ? "ensembl_gene_id" : "final_gene_name"
        cmd =
        """
        #!/usr/bin/env Rscript
        library(data.table)

        ensembl <- fread('$ensembl_file', data.table=FALSE)

        ensembl <- ensembl[ensembl\$gene_biotype %in% "protein_coding", ]
        ensembl <- ensembl[ensembl\$chromosome_name %in% c(as.character(1:22), "X"), ]

        coords <- ensembl[, c("${gene_col}", "chromosome_name", "start_position", "end_position")]
        colnames(coords) <- c("GENE", "CHR", "START", "END")
        coords <- coords[!duplicated(coords[, 1]), ]

        write.table(coords, file="gene_coords.txt", row.names=FALSE, sep="\t", col.names=TRUE, quote=FALSE)
        cat(sprintf("Written %d genes to gene_coords.txt\n", nrow(coords)))
        """
        cmd
}


// Harmonise summary statistics to LDSC format using munge_sumstats.py.
// MHC region is implicitly handled by merging against HapMap3 SNPs (w_hm3.snplist)
// and using --not-M-5-50 during the regression step.
process munge_sumstats {
    label params.ldsc.label

    container params.ldsc.container
    conda params.ldsc.conda

    publishDir "$params.rn_publish_dir/reference/ldsc/munged", mode: 'symlink'

    input:
        tuple val(pheno_id), path(sumstats_file),
              val(snp_col), val(a1_col), val(a2_col),
              val(sign_col), val(sign_null),
              val(p_col), val(n_val), val(n_case), val(n_control)
        val(reference_dir)

    output:
        tuple val(pheno_id), path("${pheno_id}.sumstats.gz"), emit: munged
        path("${pheno_id}.log", optional: true,         emit: log)

    script:
        def decompress = sumstats_file.name.endsWith(".gz")
            ? "gunzip -c ${sumstats_file} > sumstat_in.tsv"
            : "ln -s ${sumstats_file} sumstat_in.tsv"
        cmd =
        """
        ${decompress}

        munge_sumstats.py \\
            --sumstats sumstat_in.tsv \\
            --out intermediate \\
            --merge-alleles ${params.ldsc.w_hm3_snplist ?: "${reference_dir}/hapmap3/w_hm3.snplist"} \\
            --snp  ${snp_col} \\
            --a1   ${a1_col} \\
            --a2   ${a2_col} \\
            --p    ${p_col}\\
        """

        // If the sumstats are signed
        if (sign_col != "NA" && sign_null != "NA") {
            cmd += " --signed-sumstats ${sign_col},${sign_null}"
        }

        // Include the samplesize
        if (n_case != "NA" && n_control != "NA") {
            cmd += " --N-cas ${n_case} --N-con ${n_control}"
        } else if (n_val != "NA") {
            cmd += " --N ${n_val}"
        }
        
        // This was not the issue
        // Remove rows with missing values in critical columns
        cmd +=
        """
        zcat intermediate.sumstats.gz | awk -F'\t' '\$1 != "" && \$2 != "" && \$3 != "" && \$4 != "" && \$5 != ""' | gzip -c > ${pheno_id}.sumstats.gz
        """

        cmd
}


// Create per-chromosome annotation files from a gene list (make_annot.py) and
// compute LD scores for the annotation (ldsc.py --l2).
// Outputs a directory named after the condition containing all per-chromosome
// .annot.gz, .l2.ldscore.gz, .l2.M and .l2.M_5_50 files needed by ldsc_stratified.
//
// Expected reference_dir layout:
//   gene_coords.txt          – tab-separated: GENE  CHR  START  END
//   plink/{plink_prefix}{1..22}.{bed,bim,fam}
//   hapmap3/hm3_no_MHC.list   – SNP IDs to restrict LD score computation
process make_annotation_and_ldscores {
    label params.ldsc.label_high

    container params.ldsc.container
    conda params.ldsc.conda

    publishDir "$params.rn_publish_dir/ldsc/${params.rn_runname}/annotations/", mode: 'symlink'

    input:
        tuple val(condition_name), path(gene_set)
        val(reference_dir)
        path(gene_coords)

    output:
        tuple val(condition_name), path("${condition_name}/"), emit: ldscores
        path("${condition_name}_annot.log", optional: true, emit: log)

    script:
        def plink_dir   = params.ldsc.plink_dir  ?: "${reference_dir}/plink/${params.ldsc.plink_prefix}"
        def hm3_no_mhc  = params.ldsc.hm3_no_mhc ?: "${reference_dir}/hapmap3/hm3_no_MHC.list"
        cmd =
        """
        mkdir -p ${condition_name}

        for chr in \$(seq 1 22); do
            make_annot.py \\
                --gene-set-file ${gene_set} \\
                --gene-coord-file ${gene_coords} \\
                --windowsize ${params.ldsc.window_size} \\
                --bimfile ${plink_dir}\${chr}.bim \\
                --annot-file ${condition_name}/${condition_name}.\${chr}.annot.gz
        done

        for chr in \$(seq 1 22); do
            ldsc.py \\
                --l2 \\
                --bfile ${plink_dir}\${chr} \\
                --ld-wind-cm ${params.ldsc.ld_wind_cm} \\
                --annot ${condition_name}/${condition_name}.\${chr}.annot.gz \\
                --thin-annot \\
                --out ${condition_name}/${condition_name}.\${chr} \\
                --print-snps ${hm3_no_mhc}
        done

        cat .command.log > ${condition_name}_annot.log
        """
        cmd
}


// Per-chromosome variant of make_annotation_and_ldscores.
// Runs one job per chromosome so each chromosome can be parallelised across cluster nodes.
// Outputs flat files named ${condition_name}.${chr}.* (no subdirectory) so they can be
// collected with groupTuple and staged together in the ldsc_stratified working directory.
//
// Typical subworkflow usage:
//   ldscores_ch = make_annotation_and_ldscores_per_chr(
//       gene_sets_ch.combine(Channel.from(1..22)), reference_dir, gene_coords
//   ).ldscores
//       .groupTuple()                          // group all chr files by condition_name
//       .map { cond, chrs, files -> tuple(cond, files.flatten()) }
process make_annotation_and_ldscores_per_chr {
    label params.ldsc.label_high

    container params.ldsc.container
    conda params.ldsc.conda

    //publishDir "$params.rn_publish_dir/ldsc/${params.rn_runname}/annotations/", mode: 'symlink'

    input:
        tuple val(condition_name), path(gene_set), val(chr)
        val(reference_dir)
        path(gene_coords)

    output:
        tuple val(condition_name), val(chr),
              path("${condition_name}.${chr}.*"), emit: ldscores
        path("${condition_name}_${chr}_annot.log", optional: true, emit: log)

    script:
        def plink_dir   = params.ldsc.plink_dir  ?: "${reference_dir}/plink/${params.ldsc.plink_prefix}"
        def hm3_no_mhc  = params.ldsc.hm3_no_mhc ?: "${reference_dir}/hapmap3/hm3_no_MHC.list"
        cmd =
        """
        make_annot.py \\
            --gene-set-file ${gene_set} \\
            --gene-coord-file ${gene_coords} \\
            --windowsize ${params.ldsc.window_size} \\
            --bimfile ${plink_dir}${chr}.bim \\
            --annot-file ${condition_name}.${chr}.annot.gz

        ldsc.py \\
            --l2 \\
            --bfile ${plink_dir}${chr} \\
            --ld-wind-cm ${params.ldsc.ld_wind_cm} \\
            --annot ${condition_name}.${chr}.annot.gz \\
            --thin-annot \\
            --out ${condition_name}.${chr} \\
            --print-snps ${hm3_no_mhc}

        cat .command.log > ${condition_name}_${chr}_annot.log
        """
        cmd
}


// Collect flat per-chromosome ldscore files (from make_annotation_and_ldscores_per_chr)
// into a single directory, matching the output structure of make_annotation_and_ldscores
// so the result can be passed directly to ldsc_stratified.
process collect_ldscores {
    label "tiny"

    container params.ldsc.container
    conda params.ldsc.conda

    input:
        tuple val(condition_name), path(chr_files)

    output:
        tuple val(condition_name), path("${condition_name}/"), emit: ldscores

    script:
        """
        mkdir -p ${condition_name}
        cp -P ${chr_files} ${condition_name}/
        """
}


// Run partitioned (stratified) LD score regression for one phenotype x condition pair.
// The custom annotation is added on top of the baselineLD_v2.2 model.
// --overlap-annot enables enrichment estimates when annotations overlap.
//
// Expected reference_dir layout:
//   baselineLD_v2.2/baselineLD.{1..22}.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}
//   weights/weights.hm3_noMHC.{1..22}.l2.ldscore.gz
//   frq/{frq_prefix}{1..22}.frq
process ldsc_stratified {
    label params.ldsc.label_run

    container params.ldsc.container
    conda params.ldsc.conda

    publishDir "$params.rn_publish_dir/ldsc/${params.rn_runname}/results/${pheno_id}", mode: 'symlink'

    input:
        tuple val(pheno_id), path(munged_sumstats),
              val(condition_name), path(condition_ldscores)
        val(reference_dir)

    output:
        tuple val(pheno_id), val(condition_name),
              path("${pheno_id}__${condition_name}.results"), emit: results
        path("${pheno_id}__${condition_name}.log", optional: true, emit: log)

    script:
        def baseline_ld = params.ldsc.baseline_ld_chr ?: "${reference_dir}/baselineLD_v2.2/baselineLD."
        def weights     = params.ldsc.weights_chr      ?: "${reference_dir}/weights/weights.hm3_noMHC."
        def frq         = params.ldsc.frq_chr          ?: "${reference_dir}/frq/${params.ldsc.frq_prefix}"
        cmd =
        """
        ldsc.py \\
            --h2 ${munged_sumstats} \\
            --ref-ld-chr ${baseline_ld},${condition_ldscores}/${condition_name}. \\
            --w-ld-chr ${weights} \\
            --overlap-annot \\
            --thin-annot \\
            --frqfile-chr ${frq} \\
            --print-coefficients \\
            --out ${pheno_id}__${condition_name}

        cat .command.log > ${pheno_id}__${condition_name}.log
        """
        cmd
}


// Merge all per-phenotype x per-condition .results files into a single TSV.
// Columns phenotype and condition are prepended. The Prop._h2 column indicates
// what fraction of SNP heritability is captured by each gene set.
process aggregate_ldsc_results {
    label params.ldsc.label

    container params.ldsc.container
    conda params.ldsc.conda

    publishDir "$params.rn_publish_dir/ldsc/${params.rn_runname}", mode: 'symlink'

    input:
        path(results_files, arity: '1..*')

    output:
        path("ldsc_results_aggregated.tsv"), emit: aggregated
        path("ldsc_results_aggregated_annot*.tsv"), emit: per_annotation
        path("ldsc_aggregation.log", optional: true, emit: log)

    script:
        cmd =
        """
        python3 << 'PYEOF'
        import pandas as pd
        import glob
        import os
        import re

        records = []
        for fpath in sorted(glob.glob("*.results")):
            base = os.path.basename(fpath).replace(".results", "")
            parts = base.split("__", 1)
            pheno_id  = parts[0]
            condition = parts[1] if len(parts) > 1 else "unknown"
            df = pd.read_csv(fpath, sep="\\t")
            df.insert(0, "condition", condition)
            df.insert(0, "phenotype", pheno_id)
            records.append(df)

        if not records:
            raise RuntimeError("No .results files found to aggregate")

        combined = pd.concat(records, ignore_index=True)
        combined.to_csv("ldsc_results_aggregated.tsv", sep="\\t", index=False)
        print(f"Aggregated {len(records)} result files into ldsc_results_aggregated.tsv")

        annot_indices = combined["Category"].str.extract(r'_(\\d+)\$', expand=False).dropna().unique()
        for idx in sorted(annot_indices, key=int):
            subset = combined[combined["Category"].str.endswith(f"_{idx}")]
            out = f"ldsc_results_aggregated_annot{idx}.tsv"
            subset.to_csv(out, sep="\\t", index=False)
            print(f"Wrote {len(subset)} rows to {out}")
        PYEOF

        cat .command.log > ldsc_aggregation.log
        """
        cmd
}
