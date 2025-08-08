#!/usr/bin/env nextflow

process prepare_cnmf {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda

    publishDir "$params.rn_publish_dir/h5ad/cnmf", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
    output:
        path("*.Corrected.HVG.Varnorm.h5ad", emit: corrected)
        path("*.TP10K.h5ad", emit: normalized)
        path("*.Corrected.HVGs.txt", emit: hvg)
        
    script:
        cmd = "run_preprocess.py -i ${file} -o ${id}"
        
        if (params.cnmf.harmony_vars != null) {
            cmd += " --harmony_vars ${params.cnmf.harmony_vars}"
        }
        
        if (params.cnmf.seed != null) {
            cmd += " --seed ${params.cnmf.seed}"
        }
        
        if (params.cnmf.n_variable != null) {
            cmd += " --n_variable ${params.cnmf.n_variable}"
        }
        
        if (params.rn_feature_type_col != null) {
            cmd += " --feature_type_col ${params.rn_feature_type_col}"
        }
        
        cmd
}