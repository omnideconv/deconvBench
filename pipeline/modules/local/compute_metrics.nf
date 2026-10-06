// Compare the estimated cell-type fractions of one result with the ground truth of the bulk dataset.
process COMPUTE_METRICS {

      tag "${meta.id}"

      publishDir path: { "${params.results_dir_general}/${meta.result_dir}" }, mode: 'copy', saveAs: { filename -> filename.replaceFirst('^out/', '') }

      input:
      tuple val(meta), path(bulk_dir, stageAs: 'bulk/*'), path(results, stageAs: 'out/*')

      output:
      tuple val(meta), path('out/results_metric.rds'), emit: metrics

      script:
      """
      computeMetricsNF.R '${meta.sc}' 'sc' '${meta.bulk}' 'bulk' '${meta.method}' '0' '0' 'out' '${projectDir}'
      """
}
