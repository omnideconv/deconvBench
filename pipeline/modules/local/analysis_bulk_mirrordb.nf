include { publishedResult; methodSettings } from '../../functions'

// Deconvolve a real bulk dataset with each reference of the impact-technology analysis.
process ANALYSIS_BULK_MIRRORDB {

      tag "${method}_${bulk_name}_bulk_${sc_name}_signature"
      label 'process_default'

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_impact_technology}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(bulk_name), path(bulk_dir, stageAs: 'bulk/*'), val(method)

      output:
      path 'out/**', emit: results

      script:
      """
      analysisNF_mixed_simulations_real_dataset.R '${sc_name}' 'sc' '${bulk_name}' 'bulk' 'preprocess' '${method}' "\$PWD/out" '${params.ncores}' '${projectDir}'${methodSettings(method)}
      """
}
