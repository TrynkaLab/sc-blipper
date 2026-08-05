#!/usr/bin/env nextflow

// Prepare itself does not use a GPU. This variant runs on the ordinary cNMF
// queue while using the GPU-capable environment and cNMF's GPU-aware prepare
// adapter, which persists the requested solver for later GPU stages.
process cnmf_prepare_gpu {
    label params.cnmf.label
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.prepare?.container ?: params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.prepare?.conda ?: params.cnmf_gpu?.conda ?: projectDir + '/environment_gpu.yml'}"

    input:
        tuple val(id), path(file)
        path(tpm)
        path(hvg)
    output:
        tuple val(id), path("${id}")

    script:
        def seed = Math.round(Math.random() * 10000).toInteger()
        cmd =
        """
        cnmf prepare \
        --output-dir ./ \
        --name ${id} \
        -c ${file} \
        -k ${params.cnmf.k.split(",").join(' ')} \
        --n-iter ${params.cnmf.n_iter} \
        --numgenes ${params.preprocess.n_variable}\
        """

        if (params.cnmf.seed != null) {
            cmd += " --seed ${params.cnmf.seed}"
        } else {
            cmd += " --seed ${seed}"
        }

        if (tpm.getFileName().toString() != "NO_TPM") {
            cmd += " --tpm ${tpm}"
        }

        if (hvg.getFileName().toString() != "NO_HVG") {
            cmd += " --genes-file ${hvg}"
        }

        def prepare = params.cnmf_gpu?.prepare ?: [:]
        if (prepare.solver != null) {
            cmd += " --solver ${prepare.solver}"
        }

        if (prepare.beta_loss != null) {
            cmd += " --beta-loss ${prepare.beta_loss}"
        }

        cmd += " --engine gpu"

        cmd
}

process cnmf_factorize_gpu {
    label "${params.cnmf_gpu?.label ?: 'gpu_medium'}"
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.conda ?: projectDir + '/environment_gpu.yml'}"

    input:
        tuple val(id), path(file, name: "tmp/*"), val(worker_index)
        val n_workers
    output:
        //tuple val(id), path("${id}/cnmf_tmp/*.df.npz")
        path("${id}/cnmf_tmp/*.k_*.iter_*.df.npz", emit: files)

    script:
        // Process-specific GPU options override shared GPU options; unset values are left to cNMF.
        def gpu = params.cnmf_gpu ?: [:]
        def gpu_factorize = gpu.factorize ?: [:]
        def gpu_device = gpu_factorize.device != null ? gpu_factorize.device : gpu.device
        def gpu_dtype = gpu_factorize.dtype != null ? gpu_factorize.dtype : gpu.dtype
        def gpu_compile_block = gpu_factorize.compile_block != null ? gpu_factorize.compile_block : gpu.compile_block
        def gpu_allow_tf32 = gpu_factorize.allow_tf32 != null ? gpu_factorize.allow_tf32 : gpu.allow_tf32
        def gpu_compile = gpu_factorize.compile != null ? gpu_factorize.compile : gpu.compile
        def gpu_eps = gpu_factorize.eps != null ? gpu_factorize.eps : gpu.eps
        def gpu_check_every = gpu_factorize.check_every != null ? gpu_factorize.check_every : gpu.check_every
        def gpu_batch = gpu_factorize.batch != null ? gpu_factorize.batch : gpu.batch

        // Only emit explicit GPU flags so cNMF remains the source of runtime defaults.
        def gpu_device_arg = gpu_device != null ? "--gpu-device ${gpu_device}" : ""
        def gpu_dtype_arg = gpu_dtype != null ? "--gpu-dtype ${gpu_dtype}" : ""
        def gpu_compile_block_arg = gpu_compile_block != null ? "--gpu-compile-block ${gpu_compile_block}" : ""
        def gpu_allow_tf32_arg = gpu_allow_tf32 != null && gpu_allow_tf32.toString().trim().toLowerCase() in ["1", "true", "yes", "on"] ? "--gpu-allow-tf32" : ""
        def gpu_compile_arg = gpu_compile != null && gpu_compile.toString().trim().toLowerCase() in ["1", "true", "yes", "on"] ? "--gpu-compile" : ""
        def gpu_eps_arg = gpu_eps != null ? "--gpu-eps ${gpu_eps}" : ""
        def gpu_check_every_arg = gpu_check_every != null ? "--gpu-check-every ${gpu_check_every}" : ""
        def gpu_batch_arg = gpu_batch != null ? "--gpu-batch ${gpu_batch}" : ""      // factorize-only: replicates per GPU launch
        // Join only the flags that are set onto ONE line, so an unset/false optional flag can't leave a
        // blank line that breaks the shell line-continuation (previously: `--gpu-allow-tf32: command not found`).
        def gpu_extra_args = [gpu_device_arg, gpu_dtype_arg, gpu_compile_block_arg, gpu_allow_tf32_arg, gpu_compile_arg, gpu_eps_arg, gpu_check_every_arg, gpu_batch_arg].findAll { it }.join(' ')

        cmd =
        """
        mkdir -p ${id}/cnmf_tmp

        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/

        cnmf factorize \
        --output-dir ./ \
        --name ${id} \
        --worker-index ${worker_index} \
        --total-workers ${n_workers} \
        --engine gpu \
        ${gpu_extra_args}
        """

        cmd
}

process cnmf_consensus_gpu {
    label "${params.cnmf_gpu?.label ?: 'gpu_medium'}"
    scratch params.rn_scratch

    container "${params.cnmf_gpu?.container ?: params.rn_container}"
    conda "${params.cnmf_gpu?.conda ?: projectDir + '/environment_gpu.yml'}"
    publishDir "$params.rn_publish_dir/cnmf/consensus/${id}/k_${k}", mode: 'symlink', saveAs: { filepath ->
        def pathStr = filepath.toString()
        def newPathStr = pathStr.replaceFirst("^${id}/", "")
        return newPathStr
    }

    input:
        tuple val(id), path(file, name: "tmp/*"), val(k)
        path(merged, name: "merged/*")
        path(h5ad)
    output:
        tuple val(k), path("${id}/${id}.gene_spectra_score*"), emit: spectra_score
        tuple val(k), path("${id}/${id}.gene_spectra_tpm*"), emit: spectra_tpm
        tuple val(k), path("${id}/${id}.spectra.k_*"), emit: spectra_k
        tuple val(k), path("${id}/${id}.usages.k_*"), emit: usages_k
        tuple val(k), path("${id}/${id}.starcat_spectra.k_*"), emit: starcat_spectra_k
        tuple val(k), path("${id}/${id}.*.png"), emit: plots
        tuple val(k), path("${id}/${id}.*.h5ad"), optional: true, emit: h5ad
        tuple val(k), path("${id}/cnmf_tmp/${id}.local_density_cache.k_${k}.merged.df.npz"), path("${id}/cnmf_tmp/${id}.spectra.k_${k}.merged.df.npz"),  path("${id}/cnmf_tmp/${id}.nmf_idvrun_params.yaml"),  path("${id}/cnmf_tmp/${id}.norm_counts.h5ad"), emit: qc_input
    script:

        String local_dens = params.cnmf.local_density.toString().replace('.', '_')
        // Consensus can use different GPU options from factorize; unset values fall through to cNMF.
        def gpu = params.cnmf_gpu ?: [:]
        def gpu_consensus = gpu.consensus ?: [:]
        def gpu_device = gpu_consensus.device != null ? gpu_consensus.device : gpu.device
        def gpu_dtype = gpu_consensus.dtype != null ? gpu_consensus.dtype : gpu.dtype
        def gpu_compile_block = gpu_consensus.compile_block != null ? gpu_consensus.compile_block : gpu.compile_block
        def gpu_allow_tf32 = gpu_consensus.allow_tf32 != null ? gpu_consensus.allow_tf32 : gpu.allow_tf32
        def gpu_compile = gpu_consensus.compile != null ? gpu_consensus.compile : gpu.compile
        def gpu_eps = gpu_consensus.eps != null ? gpu_consensus.eps : gpu.eps
        def gpu_check_every = gpu_consensus.check_every != null ? gpu_consensus.check_every : gpu.check_every

        // Only emit explicit GPU flags so cNMF remains the source of runtime defaults.
        def gpu_device_arg = gpu_device != null ? "--gpu-device ${gpu_device}" : ""
        def gpu_dtype_arg = gpu_dtype != null ? "--gpu-dtype ${gpu_dtype}" : ""
        def gpu_compile_block_arg = gpu_compile_block != null ? "--gpu-compile-block ${gpu_compile_block}" : ""
        def gpu_allow_tf32_arg = gpu_allow_tf32 != null && gpu_allow_tf32.toString().trim().toLowerCase() in ["1", "true", "yes", "on"] ? "--gpu-allow-tf32" : ""
        def gpu_compile_arg = gpu_compile != null && gpu_compile.toString().trim().toLowerCase() in ["1", "true", "yes", "on"] ? "--gpu-compile" : ""
        def gpu_eps_arg = gpu_eps != null ? "--gpu-eps ${gpu_eps}" : ""
        def gpu_check_every_arg = gpu_check_every != null ? "--gpu-check-every ${gpu_check_every}" : ""
        // Join only the flags that are set onto ONE line, so an unset/false optional flag can't leave a
        // blank line that breaks the shell line-continuation (previously: `--gpu-allow-tf32: command not found`).
        def gpu_extra_args = [gpu_device_arg, gpu_dtype_arg, gpu_compile_block_arg, gpu_allow_tf32_arg, gpu_compile_arg, gpu_eps_arg, gpu_check_every_arg].findAll { it }.join(' ')

        cmd =
        """
        mkdir -p ${id}/cnmf_tmp

        ln -s \$(pwd)/tmp/${id}/cnmf_tmp/* ${id}/cnmf_tmp/
        ln -s \$(pwd)/tmp/${id}/*.overdispersed_genes.txt ${id}/
        ln -s \$(pwd)/merged/* ${id}/cnmf_tmp/

        cnmf consensus \
        --output-dir ./ \
        --name ${id} \
        --components ${k} \
        --local-density-threshold ${params.cnmf.local_density} \
        --show-clustering \
        --engine gpu \
        ${gpu_extra_args}

        find ${id} -maxdepth 1 -type f -name "*.txt" ! -xtype l -exec gzip {} \\;
        """

        if (h5ad.getFileName().toString() != "NO_H5AD") {
            cmd +=
            """
            cnmf_to_h5ad.py \
            --spectra ${id}/${id}.gene_spectra_score.*.txt.gz \
            --usage ${id}/${id}.usages.k_*.txt.gz \
            --obs ${h5ad} \
            --output ${id}/${id}.k_${k}.dt_${local_dens}.h5ad
            """
        }

        cmd
}
