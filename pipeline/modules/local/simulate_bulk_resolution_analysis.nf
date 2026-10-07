// Simulate pseudo-bulks with ground truth at fine, normal and (if annotated) coarse cell-type resolution.
process SIMULATE_BULK_RESOLUTION_ANALYSIS {

      tag "${sc_name}"
      label 'process_default'

      beforeScript 'chmod o+rw .'

      // simulated datasets are kept in the preprocessing directory, as before
      publishDir "${params.preProcess_dir}", mode: 'copy'

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(n_cells), val(n_samples), val(cell_types_fine), val(replicates)

      output:
      tuple val(sc_name), val(replicates), path("pseudo_bulk_resolution/${sc_name}_resolution_analysis_sim"), emit: bulk

      script:
      """
      simulateBulkNF_impact_cell_resolution.R '${sc_name}' 'sc/' '${n_cells}' '${n_samples}' '${cell_types_fine}' '${replicates}' '.' '${params.ncores}'
      """
}
