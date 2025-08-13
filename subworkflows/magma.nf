#!/usr/bin/env nextflow
include { ensembl_to_magma_geneloc; perpare_sumstats; magma_annotate; magma_prepare; magma_merge; magma_concat; magma_write_manifest } from "../processes/magma.nf"

workflow magma_base {

    take:
        magma
        convert
        rn_ensembl_version
        ensembl_reference
    main:
        // Check if recomputing magma is needed
        if (magma.manifest_magma == null) {
            
            // Read manifest
            if (magma.manifest_sumstats != null) {
                manifest_sumstats_file = file(magma.manifest_sumstats)
                if (!manifest_sumstats_file.exists()) {
                    throw new Exception("Supplied manifest_sumstats file does not exist")
                }
            } else {
                throw new Exception("Sumstat manifest is null, this is needed for magma")
            }
            
            manifest_sumstats = Channel.fromPath(magma.manifest_sumstats)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.name,
            row.n,
            row.variant_col,
            row.p_col,
            file(row.path))}
            
            // Plink reference files
            prefix = file(magma.ld_reference).name
            bedfile = file(magma.ld_reference + ".bed")
            bimfile = file(magma.ld_reference + ".bim")
            famfile = file(magma.ld_reference + ".fam")

            if (!bedfile.exists() || !bimfile.exists() || !famfile.exists()) {throw new Exception("ld_reference .bed/.bim/.fam files do not exist. Check the path, and check if its plink, and you specified only the prefix")}
                        
            // Prepare geneloc file based on ensembl with gene names or gene id's matching
            magma_geneloc = ensembl_to_magma_geneloc(ensembl_reference, !convert.invert_linker)
            
            // Prepare the magma snp to gene mapping for the reference file
            magma_gene_annot = magma_annotate("v"+rn_ensembl_version+"_ensembl", 
                                            bimfile,
                                            magma_geneloc)
            
            // Prepare the summary stats snp_pval format
            magma_sumstats = perpare_sumstats(manifest_sumstats)
            
            // ---------------------------------------------------------------
            // Prepare the magma run 
            // Make a channel per batch
            magma_prepare_in = magma_sumstats
            .map { v -> (1..magma.n_batch).collect{ i -> tuple(*v, i) } }
            .flatMap()
                        
            // Run magma prepare
            magma_prepare_out = magma_prepare(magma_prepare_in, tuple(prefix, bedfile, bimfile, famfile), magma_gene_annot)
                        
            // Groups the output and releases it as soon as it is done
            magma_merge_in = magma_prepare_out.groupTuple(by:0, size: magma.n_batch).map{row -> tuple(row[0], row[1].flatten())}
            
            // Merge the magma results over the batches
            magma_merge_out = magma_merge(magma_merge_in)
            
            // Save the output to a file for use in later runs                
            magma_write_manifest(magma_merge_out.raw.map{row -> row.join("\t") }.collect())
            
        } else {
            // Using previous magma results from manifest
            magma_merge_out = Channel.fromPath(magma.manifest_magma)
                .splitCsv(header:true, sep:"\t")
                .map { row -> tuple(
                row.name,
                row.raw)}
        }
    emit:
        raw=magma_merge_out.raw
        logs=magma_merge_out.logs
        out=magma_merge_out.out
}
    