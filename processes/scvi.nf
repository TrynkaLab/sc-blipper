#!/usr/bin/env nextflow

process scvi_cnmf {
    label params.preprocess.scvi.label
    container params.preprocess.scvi.container
    conda params.preprocess.scvi.conda

    publishDir "$params.rn_publish_dir/h5ad/scvi/${params.rn_runname}", mode: 'symlink'

    input:
        tuple val(id), path(h5ad_file)

    output:
        path("*.Corrected.HVG.Varnorm.h5ad", emit: corrected)
        path("*.TP10K.h5ad", emit: normalized)
        path("*.Corrected.HVGs.txt", emit: hvg)
        path("${id}*.png"), emit: plots, optional: true
        path("${id}*_meta_with_umap_coords.tsv"), emit: umap, optional: true
        path("${id}__scvi__nl${params.scvi.n_latent}.h5ad", emit: scvi, optional: true)
        tuple val(id), path("*.Corrected.HVG.Varnorm.h5ad"), emit: preprocessed

    script:
        cmd = 
        """
        run_scvi.py \
        --input ${h5ad_file} \
        --output_prefix ${id} \
        --batch_key ${params.scvi.batch_key} \
        --n_latent ${params.scvi.n_latent} \
        --epochs ${params.scvi.epochs} \
        --variable_genes ${params.preprocess.n_variable} \
        """

        if (params.scvi.skip_early_stopping) {
            cmd += " --skip_early_stopping"   
        }
        
        if (params.scvi.cat_covariates) {
            cmd += " --cat_covariates ${params.scvi.cat_covariates.join(' ')}"
        }
        
        if (params.scvi.cont_covariates) {
            cmd += " --cont_covariates ${params.scvi.cont_covariates.join(' ')}"
        }
        
        if (params.scvi.model_path) {
            cmd += " --model ${params.scvi.model_path}"   
        }
        
        // These are the cNMF scaling parameters
        if (params.scvi.max_scaled_thresh) {
            cmd += " --max_scaled_thresh ${params.scvi.max_scaled_thresh}"   
        }
        
        if (params.scvi.quantile_thresh) {
            cmd += " --quantile_thresh ${params.scvi.quantile_thresh}"   
        }
        
        if (params.scvi.skip_scvi_h5ad) {
            cmd += " --skip_scvi_h5ad"   
        }
        cmd
}


process scvi {
    label params.preprocess.scvi.label
    container params.preprocess.scvi.container
    conda params.preprocess.scvi.conda

    publishDir "$params.rn_publish_dir/h5ad/scvi/${params.rn_runname}", mode: 'symlink'

    input:
        tuple val(id), path(h5ad_file)

    output:
        tuple val(id), path("${id}__scvi__nl${params.scvi.n_latent}.h5ad"), emit: preprocessed
        path("${id}*.png"), emit: plots, optional: true
        path("${id}*_meta_with_umap_coords.tsv"), emit: umap, optional: true

    script:
        cmd = 
        """
        run_scvi.py \
        --input ${h5ad_file} \
        --output_prefix ${id} \
        --batch_key ${params.scvi.batch_key} \
        --n_latent ${params.scvi.n_latent} \
        --epochs ${params.scvi.epochs} \
        --variable_genes ${params.preprocess.n_variable} \
        --skip_cnmf_h5ad \
        """
                
        if (params.scvi.skip_early_stopping) {
            cmd += " --skip_early_stopping"   
        }
        
        if (params.scvi.cat_covariates) {
            cmd += " --cat_covariates ${params.scvi.cat_covariates.join(' ')}"
        }
        
        if (params.scvi.cont_covariates) {
            cmd += " --cont_covariates ${params.scvi.cont_covariates.join(' ')}"
        }
        
        if (params.scvi.model_path) {
            cmd += " --model ${params.scvi.model_path}"   
        }
        
        // These are the cNMF scaling parameters
        if (params.scvi.max_scaled_thresh) {
            cmd += " --max_scaled_thresh ${params.scvi.max_scaled_thresh}"   
        }
        
        if (params.scvi.quantile_thresh) {
            cmd += " --quantile_thresh ${params.scvi.quantile_thresh}"   
        }
        
        if (params.scvi.skip_scvi_h5ad) {
            cmd += " --skip_scvi_h5ad"   
        }
        cmd
}