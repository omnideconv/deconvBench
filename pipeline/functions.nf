// Helper functions shared by the workflows in main.nf.

// rows of optimal_normalizations.csv for the methods selected in params.method_list
def selectedMethods() {
    def table = file("${projectDir}/optimal_normalizations.csv").splitCsv(header: true)
    def known = table.collect { row -> row.method }
    def unknown = params.method_list.findAll { method -> !(method in known) }
    if (unknown) {
        error "No normalization defined in optimal_normalizations.csv for: ${unknown.join(', ')}"
    }
    return table.findAll { row -> row.method in params.method_list }
}

// directory of a dataset that is stored as <base_dir>/<name>
def datasetDir(base_dir, name) {
    return file("${base_dir}/${name}", type: 'dir', checkIfExists: true)
}

// Describe one combination of method, single-cell reference and bulk dataset.
// result_dir is the folder name that the figure scripts in visualisation/ parse.
def deconvolutionMeta(method, sc_name, bulk_name, preprocessed) {
    def name = "${method.method}_${sc_name}_${method.sc_norm}_${bulk_name}_${method.bulk_norm}"
    return [
        id: name,
        method: method.method,
        sc: sc_name,
        bulk: bulk_name,
        preprocessed: preprocessed ? 'true' : 'false',
        result_dir: preprocessed ? name : "${name}_ct0_rep0"
    ]
}

// [name, directory] of datasets that are stored as <base_dir>/<name>
def datasetChannel(names, base_dir) {
    return channel.fromList(names).map { name -> [name, datasetDir(base_dir, name)] }
}

// input of SIMULATE_BULK: every combination of simulation dataset, number of cells and number of samples
def simulationInput(scenarios, bias_types) {
    return datasetChannel(params.simulation_sc_dataset, params.data_dir_sc)
        .combine(channel.fromList(params.simulation_n_cells))
        .combine(channel.fromList(params.simulation_n_samples))
        .combine(channel.fromList(scenarios))
        .combine(channel.fromList(bias_types))
}

// input of CREATE_SIGNATURE: every combination of single-cell reference, bulk dataset and selected method
// ch_sc and ch_bulk emit [name, directory]
def deconvolutionInput(ch_sc, ch_bulk, preprocessed) {
    return ch_sc
        .combine(ch_bulk)
        .combine(channel.fromList(selectedMethods()))
        .map { sc_name, sc_dir, bulk_name, bulk_dir, method ->
            [deconvolutionMeta(method, sc_name, bulk_name, preprocessed), sc_dir, bulk_dir]
        }
}

// input of COMPUTE_METRICS: bulk dataset, runtime of the signature step and the deconvolution files
def metricsInput(ch_input, ch_signature, ch_deconvolution) {
    return ch_input
        .map { meta, _sc_dir, bulk_dir -> [meta, bulk_dir] }
        .join(ch_signature.map { meta, files -> [meta, files.find { f -> f.name == 'runtime_signature.rds' }] })
        .join(ch_deconvolution)
        .map { meta, bulk_dir, runtime_signature, deconvolution -> [meta, bulk_dir, [runtime_signature] + deconvolution] }
}

// Name under which a file from the out/ directory of a task is published.
// Scratch files of the methods (tmp/ directories) are not published.
def publishedResult(filename) {
    return filename.contains('/tmp/') ? null : filename.replaceFirst('^out/', '')
}
