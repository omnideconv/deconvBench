include { publishedResult; methodSettings } from '../../functions'

// Deconvolve a bulk dataset with a reference from which one cell type was removed.
process ANALYSIS_BULK_MISSING_CELL_TYPES {

      tag "${method}_${sc_name}_${bulk_name}_no_${cell_type}"
      label 'process_default'

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_missing_cell_types}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(bulk_name), path(bulk_dir, stageAs: 'bulk/*'), val(method), val(cell_type)

      output:
      path 'out/**', emit: results

      script:
      """
      analysisNF_missing_cell_type.R '${sc_name}' 'sc' '${bulk_name}' 'bulk' '${method}' '${cell_type}' "\$PWD/out" '${params.ncores}' '${projectDir}'${methodSettings(method)}
      """
}
