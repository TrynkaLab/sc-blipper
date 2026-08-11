class CnmfGpuConfigurer {
    static Map process(def config) {
        boolean enabled = config?.enabled == true

        [
            use_prepare  : enabled && config?.prepare?.enabled == true,
            use_factorize: enabled && config?.factorize?.enabled != false,
            use_consensus: enabled && config?.consensus?.enabled != false,
            prepare      : process_prepare_args(config),
            factorize    : process_factorize_args(config),
            consensus    : process_consensus_args(config)
        ]
    }

    private static String process_prepare_args(def config) {
        def prepare = config?.prepare

        [
            option('--solver', prepare?.solver),
            option('--beta-loss', prepare?.beta_loss),
            option('--max-nmf-iter', prepare?.max_iter)
        ].findAll { it }.join(' ')
    }

    private static String process_factorize_args(def config) {
        process_gpu_args(config, config?.factorize, true)
    }

    private static String process_consensus_args(def config) {
        process_gpu_args(config, config?.consensus, false)
    }

    private static String process_gpu_args(def shared, def process, boolean includeBatch) {
        def device = process?.device != null ? process.device : shared?.device
        def dtype = process?.dtype != null ? process.dtype : shared?.dtype
        def compileBlock = process?.compile_block != null ? process.compile_block : shared?.compile_block
        def allowTf32 = process?.allow_tf32 != null ? process.allow_tf32 : shared?.allow_tf32
        def compile = process?.compile != null ? process.compile : shared?.compile
        def eps = process?.eps != null ? process.eps : shared?.eps
        def checkEvery = process?.check_every != null ? process.check_every : shared?.check_every
        def batch = process?.batch != null ? process.batch : shared?.batch

        [
            option('--gpu-device', device),
            option('--gpu-dtype', dtype),
            option('--gpu-compile-block', compileBlock),
            booleanOption('--gpu-allow-tf32', allowTf32),
            booleanOption('--gpu-compile', compile),
            option('--gpu-eps', eps),
            option('--gpu-check-every', checkEvery),
            includeBatch ? option('--gpu-batch', batch) : ''
        ].findAll { it }.join(' ')
    }

    private static String option(String name, def value) {
        value != null ? "${name} ${value}" : ''
    }

    private static String booleanOption(String name, def value) {
        value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on'] ? name : ''
    }
}
