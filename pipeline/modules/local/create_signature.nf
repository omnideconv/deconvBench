// Build the signature (or model) of one method from one single-cell reference.
// bulk_dir is empty and meta.bulk is '-' if the signature does not depend on the bulk dataset.
// The files are published by DECONVOLUTE together with the deconvolution result.
// The script line keeps its argument order and indentation: visualisation/supplement/fig_s4.Rmd
// reads method and dataset names from the script column of the trace file by position.
process CREATE_SIGNATURE {

      tag "${meta.id}"
      label 'process_high'

      beforeScript 'chmod o+rw .'

      input:
      tuple val(meta), path(sc_dir, stageAs: 'sc/*'), path(bulk_dir, stageAs: 'bulk/*')

      output:
      tuple val(meta), path('out/*'), emit: signature

      script:
      """
      computeSignaturesNF.R '${meta.sc}' 'sc' '${meta.bulk}' 'bulk' '${meta.method}' 'out' '${meta.preprocessed}' '0' '0' '${params.ncores}' '${projectDir}'
      """
}
