class CnmfGpuConfigurer {
    static Map process(def config) {
        boolean enabled = config?.enabled == true

        [
            use_prepare  : enabled && config?.prepare?.enabled == true,
            use_factorize: enabled && config?.factorize?.enabled != false,
            use_consensus: enabled && config?.consensus?.enabled != false
        ]
    }
}
