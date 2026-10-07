include { publishedResult } from '../../functions'

// Deconvolve the pseudo-bulks of one dataset with the reference of another (or the same) dataset.
process ANALYSIS_PSEUDOBULK_MIRRORDB {

      tag "${method}_${bulk_name}_pseudobulk_${sc_name}_signature"
      label 'process_default'

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_impact_technology}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(bulk_name), val(replicates), path(bulk_dir, stageAs: 'bulk/*'), val(method)
      path(cell_types, stageAs: 'preprocess/*')

      output:
      path 'out/**', emit: results

      script:
      """
      analysisNF_mixed_simulations.R '${sc_name}' 'sc' '${bulk_name}' 'bulk/${bulk_name}' 'preprocess' '${replicates}' '${method}' "\$PWD/out" '${params.ncores}' '${projectDir}'
      """
}
