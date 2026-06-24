// Test-only smoke module: invoke cnmf's CLI path (prepare -> factorize -> combine -> consensus).
// It uses params.rn_runname as cNMF's --name, matching sc-blipper's cnmf process contract.
nextflow.enable.dsl = 2

process CNMF_SMOKE {
    conda "${projectDir}/environment.yml"

    output:
    path "${params.rn_runname}/cnmf_tmp/*.spectra.k_*.dt_*.consensus.df.npz", emit: spectra
    path "${params.rn_runname}/cnmf_tmp/*.usages.k_*.dt_*.consensus.df.npz",  emit: usages

    script:
    """
    python ${projectDir}/tests/nextflow/scripts/make_tiny_counts.py counts.h5ad

    cnmf prepare \
        --output-dir ./ \
        --name ${params.rn_runname} \
        -c counts.h5ad \
        -k ${params.cnmf_smoke.k} \
        --n-iter ${params.cnmf_smoke.n_iter} \
        --numgenes ${params.cnmf_smoke.numgenes} \
        --seed ${params.cnmf_smoke.seed}

    cnmf factorize \
        --output-dir ./ \
        --name ${params.rn_runname} \
        --worker-index ${params.cnmf_smoke.worker_index} \
        --total-workers ${params.cnmf_smoke.total_workers}

    cnmf combine \
        --output-dir ./ \
        --name ${params.rn_runname}

    cnmf consensus \
        --output-dir ./ \
        --name ${params.rn_runname} \
        --components ${params.cnmf_smoke.k} \
        --local-density-threshold ${params.cnmf_smoke.local_density_threshold}
    """
}
