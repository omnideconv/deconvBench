// Simulate one pseudo-bulk dataset from a single-cell dataset with SimBu.
process SIMULATE_BULK {

      tag "${sc_name}-ncells${n_cells}-nsamples${n_samples}-${scenario}-${bias_type}"
      label 'process_default'

      // simulated datasets are kept next to the real bulk datasets, as before
      publishDir "${params.data_dir_bulk}", mode: 'copy'

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(n_cells), val(n_samples), val(scenario), val(bias_type)

      output:
      tuple val(bulk_name), path(bulk_name), emit: bulk

      script:
      bulk_name = "${sc_name}-ncells${n_cells}-nsamples${n_samples}-${scenario}-${bias_type}-simulation"
      """
      simulateBulkNF.R '${sc_name}' 'sc' '${n_cells}' '${n_samples}' '${scenario}' '${bias_type}' '.' '${params.ncores}'
      """
}
