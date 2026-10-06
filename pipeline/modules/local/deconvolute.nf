// Deconvolve one bulk dataset with one method. The files of the signature step are staged into
// out/, the directory the script reads the signature from and writes its results to.
// The script line keeps its argument order and indentation (see create_signature.nf).
process DECONVOLUTE {

      tag "${meta.id}"
      label 'process_default'

      publishDir path: { "${params.results_dir_general}/${meta.result_dir}" }, mode: 'copy', saveAs: { filename -> filename.replaceFirst('^out/', '') }

      input:
      tuple val(meta), path(sc_dir, stageAs: 'sc/*'), path(bulk_dir, stageAs: 'bulk/*'), path(signature, stageAs: 'out/*')

      output:
      tuple val(meta), path('out/{deconvolution.rds,runtime_deconvolution.rds}'), emit: deconvolution

      script:
      """
	runDeconvolutionNF.R '${meta.sc}' 'sc' '${meta.bulk}' 'bulk' '${meta.method}' 'out' '${meta.preprocessed}' '0' '0' '${params.species_sc}' '${params.ncores}' '${projectDir}'
	"""
}
