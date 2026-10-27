// Simulate pseudo-bulks that mirror the cell-type composition of each single-cell dataset,
// restricted to the cell types that all datasets have in common.
// The script simulates every dataset in sc_dirs; one task is run per dataset and keeps that
// dataset's pseudo-bulks.
process SIMULATE_PSEUDOBULK_MIRRORDB {

      tag "${sc_name}"
      label 'process_high'

      beforeScript 'chmod o+rw .'

      // simulated datasets are kept in the preprocessing directory, as before
      publishDir "${params.preProcess_dir}", mode: 'copy'

      input:
      tuple val(sc_name), val(n_cells), val(n_samples), val(replicates)
      path(sc_dirs, stageAs: 'sc/*')

      output:
      tuple val(sc_name), val(replicates), path("pseudo_bulk_impact_technology/${sc_name}"), emit: bulk
      path 'cell_types.rds', emit: cell_types

      script:
      """
      simulateBulkNF_mirrordb.R '${sc_name}' 'sc' '${n_cells}' '${n_samples}' '${params.sc_datasets_impact_technology}' '${replicates}' '.' '${params.ncores}'
      """
}
