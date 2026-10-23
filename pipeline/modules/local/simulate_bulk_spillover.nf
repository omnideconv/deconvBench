// Simulate pseudo-bulks that each contain a single cell type.
process SIMULATE_BULK_SPILLOVER {

      tag "${sc_name}"
      label 'process_default'

      beforeScript 'chmod o+rw .'

      // simulated datasets are kept in the preprocessing directory, as before
      publishDir "${params.preProcess_dir}", mode: 'copy'

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(n_cells), val(n_samples)

      output:
      tuple val("${sc_name}_spillover_sim"), path("pseudo_bulk_spillover/${sc_name}_spillover_sim"), emit: bulk

      script:
      """
      simulateBulkNF_spillover_analysis.R '${sc_name}' 'sc' '${n_cells}' '${n_samples}' '${params.spillover_celltypes}' '.' '${params.ncores}'
      """
}
