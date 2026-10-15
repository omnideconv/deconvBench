// Draw a subsample of a single-cell dataset (subsampling workflow).
process PREPROCESS_SINGLE_CELL {

      tag "${sc_name}_perc${ct_fraction}_rep${replicate}"
      label 'process_default'

      // subsampled datasets are kept in the preprocessing directory, as before
      publishDir "${params.preProcess_dir}", mode: 'copy'

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*'), val(ct_fraction), val(replicate)

      output:
      path "${sc_name}_perc*_rep${replicate}", type: 'dir', emit: subsample

      script:
      """
      preprocessSingleCellNF.R '${sc_name}' 'sc' '${ct_fraction}' '${replicate}' '.' '${projectDir}'
      """
}
