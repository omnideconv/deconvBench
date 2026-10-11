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

// File whose name carries the state of the files of a dataset (paths, sizes, modification times).
// Nextflow notices a change of a directory input only when files are added to or removed from it.
// With this file among the inputs, -resume also re-runs the tasks of a dataset in which a file
// was replaced.
def datasetStamp(dir) {
    def entries = []
    dir.eachFileRecurse { f -> entries << "${dir.relativize(f)}:${f.size()}:${f.lastModified()}".toString() }
    def listing = entries.sort().join('\n') + '\n'
    def stamp = file("${workflow.workDir}/dataset_stamps/${dir.name}.${listing.md5()}.stamp")
    if (!stamp.exists()) {
        stamp.parent.mkdirs()
        stamp.text = listing
    }
    return stamp
}

// [name, [directory, stamp file]] of datasets that are stored as <base_dir>/<name>
// The processes stage both next to each other, e.g. sc/<name>/ and sc/<name>.<hash>.stamp.
def datasetChannel(names, base_dir) {
    return channel.fromList(names).map { name ->
        def dir = datasetDir(base_dir, name)
        [name, [dir, datasetStamp(dir)]]
    }
}

// Is a dataset that is used for the simulation and as reference split by donor?
def splitEnabled() {
    return params.containsKey('split_reference') && params.split_reference
}

// value of a split parameter, or its default if the config does not set it
def splitSetting(name) {
    def defaults = [split_fraction_simulation: 0.5, split_min_cells: 50, split_min_cells_action: 'error', split_seed: 1, split_reference_cells_fraction: 1]
    return params.containsKey(name) && params[name] != null ? params[name] : defaults[name]
}

// folder of a split in <preProcess_dir>/reference_split; the same parameters give the same folder
def splitFolder(sc_name) {
    return "${sc_name}_simulation${splitSetting('split_fraction_simulation')}_mincells${splitSetting('split_min_cells')}_seed${splitSetting('split_seed')}_cells${splitSetting('split_reference_cells_fraction')}"
}

// Names of the datasets from which pseudo-bulks are simulated. With split_reference the
// references are used if simulation_sc_dataset is not set.
def simulationDatasets() {
    if (params.containsKey('simulation_sc_dataset') && params.simulation_sc_dataset != null) {
        return params.simulation_sc_dataset
    }
    if (splitEnabled() && params.containsKey('single_cell_list') && params.single_cell_list != null) {
        return params.single_cell_list
    }
    error "The config does not set: simulation_sc_dataset"
}

// names of the datasets that have both roles and are therefore split; none without split_reference
def datasetsToSplit(simulation_names, reference_names) {
    if (!splitEnabled()) {
        return []
    }
    requireParams(['preProcess_dir'])
    return simulation_names.intersect(reference_names)
}

// input of SIMULATE_BULK: every combination of simulation dataset, number of cells and number of samples
// ch_simulation emits [name, directory]
def simulationInput(ch_simulation, scenarios, bias_types) {
    return ch_simulation
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

// Settings of a method, passed to its script as last argument, e.g. " 'marker_fraction=0.1'".
// Empty for methods without settings, so that their commands stay as they are. Being part of the
// command, a changed setting re-runs the tasks of the method.
def methodSettings(method) {
    if (method == 'music-hc') {
        def fraction = params.containsKey('music_hc_marker_fraction') ? params.music_hc_marker_fraction : 0.1
        return " 'marker_fraction=${fraction}'"
    }
    return ''
}

// Name under which a file from the out/ directory of a task is published.
// Scratch files of the methods (tmp/ directories) are not published.
def publishedResult(filename) {
    return filename.contains('/tmp/') ? null : filename.replaceFirst('^out/', '')
}
