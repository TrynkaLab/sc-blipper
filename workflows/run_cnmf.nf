#!/usr/bin/env nextflow

include { seurat_to_h5ad; link_h5ad; } from '../processes/convert_to_h5ad.nf'
include { merge_h5ad } from '../processes/merging.nf'
include { fetch_gene_id_reference; invert_id_link } from "../processes/utils.nf"
include { cnmf_pre_process; cnmf_prepare; cnmf_factorize; cnmf_combine; cnmf_kselection; cnmf_consensus } from "../processes/cnmf.nf"
include { gsea; ora; decoupler } from "../processes/enrichment.nf"
include { ensembl_to_magma_geneloc; perpare_sumstats; magma_annotate; magma_prepare; magma_merge; magma_assoc; magma_concat } from "../processes/magma.nf"

// Main workflow
workflow cnmf {
    //TODO: make gene linker file ids unique
    main:
        // Always fetch this even if it is not used, can make it optional later
        ensembl_reference = fetch_gene_id_reference(params.rn_ensembl_version)    
        
        // If needed, fetch ensembl file
        if (params.convert.convert_gene_names) {
            if (params.convert.id_linker != null) {
                id_linker = file(params.convert.id_linker)
                if (!id_linker.exists()) {
                    throw new Exception("Supplied id linker file does not exist")
                }
                // Custom conversion file
                id_linker = Channel.value(file(params.convert.id_linker))
                id_linker_inv = invert_id_link(id_linker)
            } else {        
                        
                if (params.convert.invert_linker) {
                    // ENSEMBL > gene name conversion
                    id_linker = ensembl_reference.ensembl_to_name
                    id_linker_inv = ensembl_reference.name_to_ensembl
                } else {
                    // gene name > ENSEMBL conversion
                    id_linker = ensembl_reference.name_to_ensembl
                    id_linker_inv = ensembl_reference.ensembl_to_name
                }
            }
        } else {
            // In case there is no mapping file
            id_linker = Channel.value(file("NO_MAPPING"))
            id_linker_inv = Channel.value(file("NO_MAPPING"))
        }
    
        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.id,
            file(row.file),
            row.convert_ids.equalsIgnoreCase("true"))}
        
        // Auto detect if it is a h5ad already or a seurat file
        manifest_split = manifest.branch{ it ->
            h5ad: it[1].extension == "h5ad"
            seu: ['rds', 'RDS', 'Rds'].contains(it[1].extension)
            other: false
        }
        
        // Converts seurat files to h5ad
        convert_out_a = seurat_to_h5ad(manifest_split.seu, id_linker).h5ad

        // Converts just symplink h5ad the files so they are in the ouput
        // In case conversion is needed, do that here
        convert_out_b = link_h5ad(manifest_split.h5ad, id_linker).h5ad

        // Merge the channels, they should now contain all h5ad files
        convert_out = convert_out_a.concat(convert_out_b)
        convert_out_flat = convert_out.flatMap{row -> {row[1]}}.collect()

        // Merge the .h5ad files into a single file
        merge_out = merge_h5ad(params.rn_runname, convert_out_flat).merged
        
        //--------------------------------------------------------
        // Optionally run cNMF preprocessing, which creates harmony corrected counts
        if (params.cnmf.preprocess) {
            // Preprocess the cnmf file
            cnmf_preprocess = cnmf_pre_process(merge_out)
            
            // Channel with run name and .h5ad file with counts
            cnmf_in = cnmf_preprocess.preproccessed
        } else {
            // Otherwise use the merged file, in case there is one batch,
            // just uses the other file
            cnmf_in = merge_out
        }
            
        //--------------------------------------------------------
        //                       cNMF
        //--------------------------------------------------------
        // Prepare the cnmf input folder 
        cnmf_prepared = cnmf_prepare(cnmf_in)
        
        // Create the channel with the n_workers for jobs to run in parallel
        cnmf_factorize_in = cnmf_prepared
         .map { v -> (0..(params.cnmf.n_workers-1)).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        // Perform the factorizations
        cnmf_factorize_out = cnmf_factorize(cnmf_factorize_in)
        
        // Wait for all the workers to complete, then:
        // Flatten all the output into a single path object
        cnmf_factorize_out = cnmf_factorize_out.collect()
        
        // Combine the output   
        cnmf_combine_out = cnmf_combine(cnmf_prepared, cnmf_factorize_out)     
 
        // Make the k selection plot
        cnmf_kselection(cnmf_prepared, cnmf_combine_out)
        
        // Make the consensus factorizations
        cnmf_consensus_in = cnmf_prepared
         .map { v -> (params.cnmf.k.split(",")).collect{ i -> tuple(*v, i) } }
         .flatMap()
        
        cnmf_out = cnmf_consensus(cnmf_consensus_in, cnmf_combine_out)
        // This is the end of the cnmf processing

        //--------------------------------------------------------
        // Optional post-proccessing of cnmf
        //--------------------------------------------------------

        //----------------------------------------------------------------------------------
        // Run gene set enrichment
        if (params.cnmf.run_enrichment && params.enrich.gmt_files != null) {
            
            // Channel with gmt files
            gmt_files = Channel.from(params.enrich.gmt_files.split(",")).collect()

            // Universe channel
            if (params.enrich.universe) {
                universe = file(params.enrich.universe)
                if (!id_linker.exists()) {
                    throw new Exception("Supplied id universe file does not exist")
                }
                // Custom conversion file
                universe = Channel.value(file(params.enrich.universe))
            } else {
                // In case there is no universe file
                universe = Channel.value(file("NO_UNIVERSE"))
            }
    
            if (params.cnmf.k_ignore != null) {
                gsea_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                gsea_in = cnmf_out.spectra_k.map{i -> tuple("k_"+i[0], i[1])}
            }

            // Run GSEA
            gsea_out = gsea("cnmf/consensus/${params.rn_runname}",
                            gsea_in,
                            gmt_files,
                            universe,
                            true)
                            
            // Run ORA
            ora_out = ora("cnmf/consensus/${params.rn_runname}",
                gsea_in,
                gmt_files,
                universe,
                true,
                false,
                0,
                params.enrich.use_top)
        }
        
        //----------------------------------------------------------------------------------
        // Run decoupler (Progeny + collecTRI)
        if (params.cnmf.run_decoupler) {
            
            if (params.cnmf.k_ignore != null) {
                decoupler_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
            } else {
                decoupler_in = cnmf_out.spectra_k.map{i -> tuple("k_"+i[0], i[1])}
            }
            
            if (params.convert.invert_linker) {
                // This assumes you have converted to gene symbols
                decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, file("NO_MAPPING"))
            } else {
                // In the case you converted everything to ensembl names, keep things consistent and convert progeny as well
                decoupler_out = decoupler("cnmf/consensus/${params.rn_runname}", decoupler_in, true, id_linker_inv) 
            }
        }
        
        //----------------------------------------------------------------------------------
        // Run magma
        if (params.cnmf.run_magma) {
            
            // Check if recomputing magma is needed
            if (params.magma.manifest_magma == null) {
                
                // Read manifest
                if (params.magma.manifest_sumstats != null) {
                    manifest_sumstats_file = file(params.magma.manifest_sumstats)
                    if (!manifest_sumstats_file.exists()) {
                        throw new Exception("Supplied manifest_sumstats file does not exist")
                    }
                } else {
                    throw new Exception("Sumstat manifest is null, this is needed for magma")
                }
                
                manifest_sumstats = Channel.fromPath(params.magma.manifest_sumstats)
                .splitCsv(header:true, sep:"\t")
                .map { row -> tuple(
                row.name,
                row.n,
                row.variant_col,
                row.p_col,
                file(row.path))}
                
                // Plink reference files
                prefix = file(params.magma.ld_reference).name
                bedfile = file(params.magma.ld_reference + ".bed")
                bimfile = file(params.magma.ld_reference + ".bim")
                famfile = file(params.magma.ld_reference + ".fam")

                if (!bedfile.exists() || !bimfile.exists() || !famfile.exists()) {throw new Exception("ld_reference .bed/.bim/.fam files do not exist. Check the path, and check if its plink, and you specified only the prefix")}
                            
                // Prepare geneloc file based on ensembl with gene names or gene id's matching
                magma_geneloc = ensembl_to_magma_geneloc(ensembl_reference.ensembl, !params.convert.invert_linker)
                
                // Prepare the magma snp to gene mapping for the reference file
                magma_gene_annot = magma_annotate("v"+params.rn_ensembl_version+"_ensembl", 
                                                bimfile,
                                                magma_geneloc)
                
                // Prepare the summary stats snp_pval format
                magma_sumstats = perpare_sumstats(manifest_sumstats)
                
                // ---------------------------------------------------------------
                // Prepare the magma run 
                // Make a channel per batch
                magma_prepare_in = magma_sumstats
                .map { v -> (1..params.magma.n_batch).collect{ i -> tuple(*v, i) } }
                .flatMap()
                            
                // Run magma prepare
                magma_prepare_out = magma_prepare(magma_prepare_in, tuple(prefix, bedfile, bimfile, famfile), magma_gene_annot)
                            
                // Groups the output and releases it as soon as it is done
                magma_merge_in = magma_prepare_out.groupTuple(by:0, size: params.magma.n_batch).map{row -> tuple(row[0], row[1].flatten())}
                
                // Merge the magma results over the batches
                magma_merge_out = magma_merge(magma_merge_in)
            } else {
                // Using previous magma results from manifest
                magma_merge_out = Channel.fromPath(params.magma.manifest_magma)
                    .splitCsv(header:true, sep:"\t")
                    .map { row -> tuple(
                    row.name,
                    row.raw)}
            }
            
            // Magma association with cnmf
            magma_cnmf_in = cnmf_out.spectra_k.filter { tuple ->
                    def (k, path) = tuple
                    !(k in params.cnmf.k_ignore.split(','))
                }.map{i -> tuple("k_"+i[0], i[1])}
                
            magma_assoc_in = magma_merge_out.raw.combine(magma_cnmf_in)
            
            // Run regression based magma
            magma_assoc_out = magma_assoc(magma_assoc_in, true)
         
            // Concat the results in a single table
            magma_concat(params.rn_runname, magma_assoc_out.out.collect())
            
        }
        
        
        // Collect enrichment results for a summary

}