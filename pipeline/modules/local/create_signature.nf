// Build the signature (or model) of one method from one single-cell reference.
// The script line keeps its argument order and indentation: visualisation/supplement/fig_s4.Rmd
// reads method and dataset names from the script column of the trace file by position.
process CREATE_SIGNATURE {

      tag "${meta.id}"
      label 'process_high'

      beforeScript 'chmod o+rw .'

      publishDir path: { "${params.results_dir_general}/${meta.result_dir}" }, mode: 'copy', saveAs: { filename -> filename.replaceFirst('^out/', '') }

      input:
      tuple val(meta), path(sc_dir, stageAs: 'sc/*'), path(bulk_dir, stageAs: 'bulk/*')

      output:
      tuple val(meta), path('out/{signature.rds,runtime_signature.rds,model}'), emit: signature
      // files that DWLS writes while building its signature
      path 'out/*.{csv,RData}', optional: true, emit: intermediates


      script:
      """
      computeSignaturesNF.R '${meta.sc}' 'sc' '${meta.bulk}' 'bulk' '${meta.method}' 'out' '${meta.preprocessed}' '0' '0' '${params.ncores}' '${projectDir}'
      """
}
