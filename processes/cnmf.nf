#!/usr/bin/env nextflow

process cnmf_pre_process {
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
        tuple val(id), path("*.Corrected.HVG.Varnorm.h5ad"), emit: preproccessed
        
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

process cnmf_prepare {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda

    publishDir "$params.rn_publish_dir/cnmf/prepare/", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
    output:
        tuple val(id), path("${id}")
        
    script:
        cmd = 
        """
        cnmf prepare \
        --output-dir ./ \
        --name $id \
        -c $file \
        -k ${params.cnmf.k.join(' ')} \
        --n-iter ${params.cnmf.n_iter} \
        --seed ${params.cnmf.seed} \
        --numgenes ${params.cnmf.n_variable}
        """
    
}

process cnmf_factorize {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    //publishDir "$params.rn_publish_dir/cnmf/${id}/factorize", mode: 'symlink'
    
    input:
        tuple val(id), path(file), val(worker_index)
    output:
        //tuple val(id), path("${id}/cnmf_tmp/*.df.npz")
        path("${id}/cnmf_tmp/*.k_*.iter_*.df.npz", emit: files)
        
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Make a tmpdir to store the symlinks to the output of the prepare process
        mkdir -p tmp
        mv ${file} tmp/
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/

        # Now we can run the process, so the output is written to the ${id}/cnmf_tmp in the workdir, 
        # not the symlinked one
        cnmf factorize \
        --output-dir ./ \
        --name ${id} \
        --worker-index ${worker_index} \
        --total-workers ${params.cnmf.n_workers}
        """
    
        cmd
}

process cnmf_combine {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    //publishDir "$params.rn_publish_dir/cnmf/${id}/combine", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
        path(factorize)
    output:
        //tuple val(id), path("${id}")
        path("${id}/cnmf_tmp/*.k_*.merged.df.npz", emit: files)
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Make a tmpdir to store the symlinks to the output of the prepare process
        mkdir -p tmp
        mv ${file} tmp/
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${file}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${file}/*.overdispersed_genes.txt ${id}/
        mv ${factorize} ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        cnmf combine \
        --output-dir ./ \
        --name ${id}
        """
    
        cmd
}

process cnmf_kselection {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    publishDir "$params.rn_publish_dir/cnmf/consensus", mode: 'symlink'
    
    input:
        tuple val(id), path(file)
        path(merged)
    output:
        path("${id}/${id}.k_selection*")
        
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Make a tmpdir to store the symlinks to the output of the prepare process
        mkdir -p tmp
        mv ${file} tmp/
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${file}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${file}/*.overdispersed_genes.txt ${id}/
        mv ${merged} ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        
        # K selection plot for cNMF
        cnmf k_selection_plot --output-dir ./ --name ${id}
        """
    
        cmd
}


process cnmf_consensus {
    label params.cnmf.label
    scratch params.rn_scratch
    
    container params.cnmf.container
    conda params.cnmf.conda
    // I dont think this needs to be publised long term, but for now its handy for debugging
    publishDir "$params.rn_publish_dir/cnmf/consensus/k_${k}/", mode: 'symlink'
    
    input:
        tuple val(id), path(file), val(k)
        path(merged)
    output:
        path("${id}/${id}*")
        
    script:
        cmd = 
        """
        # We need to do some jujitsu to get this to work, as nextflow cannot write outside the process dir
        
        # Make a tmpdir to store the symlinks to the output of the prepare process
        mkdir -p tmp
        mv ${file} tmp/
        
        # Create an output folder in the workdir with the same name
        mkdir -p ${id}/cnmf_tmp
        
        # Symlink the files, not the whole folder, this will trick cNMF into thinking it is the same folder
        ln -s \$(pwd)/tmp/${file}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${file}/*.overdispersed_genes.txt ${id}/
        mv ${merged} ${id}/cnmf_tmp/
        
        # Now we can run the process, so the output is written to the ${id} in the workdir, 
        # not the symlinked one
        
        # Consensus for cNMF
        cnmf consensus \
        --output-dir ./ \
        --name ${id} \
        --components ${k} \
        --local-density-threshold ${params.cnmf.local_density} \
        --show-clustering
        """
    
        cmd
}
