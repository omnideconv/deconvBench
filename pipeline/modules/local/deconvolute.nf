include { publishedResult; methodSettings } from '../../functions'

// Deconvolve one bulk dataset with one method and publish the result directory.
// The files of the signature step are staged in signature/ and copied into the result by the
// script. signature is empty for methods that do not build a signature.
// The script line keeps its argument order and indentation (see create_signature.nf).
process DECONVOLUTE {

      tag "${meta.id}"
      label 'process_default'

      publishDir path: { "${params.results_dir_general}/${meta.result_dir}" }, mode: 'copy', saveAs: { filename -> publishedResult(filename) }

      input:
      tuple val(meta), path(sc_dir, stageAs: 'sc/*'), path(bulk_dir, stageAs: 'bulk/*'), path(signature, stageAs: 'signature/*')

      output:
      tuple val(meta), path('out/**'), emit: results

      script:
      """
	runDeconvolutionNF.R '${meta.sc}' 'sc' '${meta.bulk}' 'bulk' '${meta.method}' 'out' '${meta.preprocessed}' '0' '0' '${params.species_sc}' '${params.ncores}' '${projectDir}' 'signature'${methodSettings(meta.method)}
	"""
}
