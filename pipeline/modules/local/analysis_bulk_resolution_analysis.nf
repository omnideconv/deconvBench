include { publishedResult } from '../../functions'

// Deconvolve the pseudo-bulks of SIMULATE_BULK_RESOLUTION_ANALYSIS with the reference annotated at each resolution.
process ANALYSIS_BULK_RESOLUTION_ANALYSIS {

      tag "${method}_${sc_name}"
      label 'process_high'

      beforeScript 'chmod o+rw .'

      publishDir "${params.results_dir_resolution}", mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(sc_name), val(replicates), path(bulk_dir, stageAs: 'bulk/*'), path(sc_dir, stageAs: 'sc/*'), val(method)

      output:
      path 'out/**', emit: results

      script:
      bulk_name = "${sc_name}_resolution_analysis_sim"
      """
      analysisNF_impact_cell_resolution.R '${sc_name}' 'sc' '${bulk_name}' 'bulk/${bulk_name}' '${method}' '${params.cell_types_finer_res}' '${replicates}' "\$PWD/out" '${params.ncores}' '${projectDir}'
      """
}
