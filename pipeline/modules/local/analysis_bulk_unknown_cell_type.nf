include { publishedResult } from '../../functions'

// Deconvolve the pseudo-bulks of SIMULATE_BULK_UNKNOWN_CELL_TYPE with a reference of the known cell types.
process ANALYSIS_BULK_UNKNOWN_CELL_TYPE {

      tag "${method}_${sc_name}_${unknown_cell_type}"

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_unknown_content}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(sc_name), val(unknown_cell_type), val(replicates), path(bulk_dir, stageAs: 'bulk/*'), path(sc_dir, stageAs: 'sc/*'), val(method)

      output:
      path 'out/**', emit: results

      script:
      bulk_name = "${sc_name}_unknown_content_sim"
      """
      analysisNF_unknown_content.R '${sc_name}' 'sc' '${bulk_name}' 'bulk/${bulk_name}' '${replicates}' '${method}' '${params.fractions_unknown_cells}' '${params.known_cell_types}' '${unknown_cell_type}' "\$PWD/out" '${params.ncores}' '${projectDir}'
      """
}
