// Helper functions shared by the workflows in main.nf.

// Stop with a clear message if the config does not set a parameter that the workflow needs.
def requireParams(names) {
    def missing = names.findAll { name -> !params.containsKey(name) || params[name] == null }
    if (missing) {
        error "The config does not set: ${missing.join(', ')}"
    }
}

// rows of optimal_normalizations.csv for the methods selected in params.method_list
// signature_input says what the signature of a method is built from:
//   none     the method builds no signature
//   sc       single-cell reference only; one signature serves all bulk datasets
//   sc+bulk  single-cell reference and bulk dataset
def selectedMethods() {
    def table = file("${projectDir}/optimal_normalizations.csv").splitCsv(header: true)
    def known = table.collect { row -> row.method }
    def unknown = params.method_list.findAll { method -> !(method in known) }
    if (unknown) {
        error "No entry in optimal_normalizations.csv for: ${unknown.join(', ')}"
    }
    def invalid = table.findAll { row -> !(row.signature_input in ['none', 'sc', 'sc+bulk']) }
    if (invalid) {
        error "signature_input in optimal_normalizations.csv has to be none, sc or sc+bulk: ${invalid.collect { row -> row.method }.join(', ')}"
    }
    return table.findAll { row -> row.method in params.method_list }
}

// directory of a dataset that is stored as <base_dir>/<name>
def datasetDir(base_dir, name) {
    return file("${base_dir}/${name}", type: 'dir', checkIfExists: true)
}

// Describe one combination of method, single-cell reference and bulk dataset.
// result_dir is the folder name that the figure scripts in visualisation/ parse.
// signature_id names the signature this combination uses (null: the method builds none). With
// share_signatures, a signature that does not depend on the bulk dataset has the same id for all
// bulk datasets and is built once.
def deconvolutionMeta(method, sc_name, bulk_name, preprocessed, share_signatures) {
    def name = "${method.method}_${sc_name}_${method.sc_norm}_${bulk_name}_${method.bulk_norm}"
    def signature_id = name
    if (share_signatures && method.signature_input == 'none') {
        signature_id = null
    } else if (share_signatures && method.signature_input == 'sc') {
        signature_id = "${method.method}_${sc_name}_${method.sc_norm}"
    }
    return [
        id: name,
        method: method.method,
        sc: sc_name,
        bulk: bulk_name,
        preprocessed: preprocessed ? 'true' : 'false',
        result_dir: preprocessed ? name : "${name}_ct0_rep0",
        signature_id: signature_id
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

// every combination of single-cell reference, bulk dataset and selected method: [meta, sc_dir, bulk_dir]
// ch_sc and ch_bulk emit [name, directory]
def deconvolutionInput(ch_sc, ch_bulk, preprocessed, share_signatures) {
    return ch_sc
        .combine(ch_bulk)
        .combine(channel.fromList(selectedMethods()))
        .map { sc_name, sc_dir, bulk_name, bulk_dir, method ->
            [deconvolutionMeta(method, sc_name, bulk_name, preprocessed, share_signatures), sc_dir, bulk_dir]
        }
}

// input of CREATE_SIGNATURE: one entry per signature that has to be built
def signatureInput(ch_input) {
    return ch_input
        .filter { meta, _sc_dir, _bulk_dir -> meta.signature_id != null }
        .map { meta, sc_dir, bulk_dir ->
            def shared = meta.signature_id != meta.id
            def signature_meta = [id: meta.signature_id, method: meta.method, sc: meta.sc, bulk: shared ? '-' : meta.bulk, preprocessed: meta.preprocessed]
            [signature_meta, sc_dir, shared ? [] : bulk_dir]
        }
        .unique { signature_meta, _sc_dir, _bulk_dir -> signature_meta.id }
}

// input of DECONVOLUTE: every combination with the files of its signature (empty if the method builds none)
// ch_signature emits [signature_meta, files]
def deconvolutionWithSignature(ch_input, ch_signature) {
    def ch_with_signature = ch_input
        .filter { meta, _sc_dir, _bulk_dir -> meta.signature_id != null }
        .map { meta, sc_dir, bulk_dir -> [meta.signature_id, meta, sc_dir, bulk_dir] }
        .combine(ch_signature.map { signature_meta, files -> [signature_meta.id, files] }, by: 0)
        .map { _signature_id, meta, sc_dir, bulk_dir, files -> [meta, sc_dir, bulk_dir, files] }
    def ch_without_signature = ch_input
        .filter { meta, _sc_dir, _bulk_dir -> meta.signature_id == null }
        .map { meta, sc_dir, bulk_dir -> [meta, sc_dir, bulk_dir, []] }
    return ch_with_signature.mix(ch_without_signature)
}

// input of COMPUTE_METRICS: bulk dataset, deconvolution result and the runtimes of both steps
// ch_results emits [meta, all files of the result directory]
def metricsInput(ch_input, ch_results) {
    def needed = ['deconvolution.rds', 'runtime_signature.rds', 'runtime_deconvolution.rds']
    return ch_input
        .map { meta, _sc_dir, bulk_dir -> [meta, bulk_dir] }
        .join(ch_results.map { meta, files -> [meta, files.findAll { f -> f.name in needed && f.parent.name == 'out' }] })
}

// Name under which a file from the out/ directory of a task is published.
// Scratch files of the methods (tmp/ directories) are not published.
def publishedResult(filename) {
    return filename.contains('/tmp/') ? null : filename.replaceFirst('^out/', '')
}
