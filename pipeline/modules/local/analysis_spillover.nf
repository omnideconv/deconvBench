include { publishedResult } from '../../functions'

// Deconvolve the single-cell-type pseudo-bulks with the full reference.
process ANALYSIS_SPILLOVER {

      tag "${method}_${sc_name}_${bulk_name}"
      label 'process_high'

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_spillover}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(bulk_name), path(bulk_dir, stageAs: 'bulk/*'), val(method), val(sc_name), path(sc_dir, stageAs: 'sc/*')

      output:
      path 'out/**', emit: results

      script:
      """
      analysisNF_spillover.R '${sc_name}' 'sc' '${bulk_name}' 'bulk/${bulk_name}' '${method}' '${params.spillover_celltypes}' "\$PWD/out" '${params.ncores}' '${projectDir}'
      """
}
