// Simulate pseudo-bulks with increasing fractions of a cell type that will be missing from the reference.
process SIMULATE_BULK_UNKNOWN_CELL_TYPE {

      tag "${sc_name}_${unknown_cell_type}"

      beforeScript 'chmod o+rw .'

      // simulated datasets are kept in the preprocessing directory, as before
      publishDir "${params.preProcess_dir}", mode: 'copy'

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(n_cells), val(n_samples), val(unknown_cell_type), val(replicates)

      output:
      tuple val(sc_name), val(unknown_cell_type), val(replicates), path("pseudo_bulk_unknown_content/${sc_name}_unknown_content_sim"), emit: bulk

      script:
      """
      simulateBulkNF_unknown_content_analysis.R '${sc_name}' 'sc' '${n_cells}' '${n_samples}' '${params.fractions_unknown_cells}' '${params.known_cell_types}' '${unknown_cell_type}' '${replicates}' '.' '${params.ncores}'
      """
}
