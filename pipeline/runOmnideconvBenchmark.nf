#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process createSignature {

	publishDir params.results_dir

	input:
	each sc
	each rnaseq
	each mode 
	path data

	output:
	tuple val(sc), val(rnaseq), val(mode), file('*.rds')

	"""
	echo '$sc $data $mode $rnaseq'
	Rscript /nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/computeSignaturesNF.R --path '$data' --single_cell '$sc' --rna_seq '$rnaseq' --method '$mode'
	""" 
}

process test { //this is only for testing using createSignatures output

	publishDir params.results_dir

	input:
	tuple val(sc), val(rnaseq), val(mode), val(signature) 
	path data
	//how can i use all tuples with this path and not only one??
	//how can i access the signatures?? when i load it as file, i get an error: file not found since it's only the filename and not the work/.../whatever path
	
	output:
	stdout

	"""
	echo '$sc $rnaseq $signature $data'
	""" 
}

process deconvolute {

	publishDir params.results_dir

	input:
	val sc
	each rnaseq
	each mode
	path data
	file signature

	output:
	path '*.rds'

	"""
	echo '$sc $data $mode $rnaseq $signature'
	Rscript /nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/runDeconvolutionNF.R --path '$data' --single_cell '$sc' --rna_seq '$rnaseq' --method '$mode' --signature '$signature'
	""" 
}


workflow{
  single_cell = Channel.fromList(params.single_cell_list)
  rna_seq = Channel.fromList(params.rnaseq_list)
  methods = params.method_list
  data = Channel.fromPath(params.data_dir)
  createSignature(single_cell, rna_seq, methods, data)
  test(createSignature.out, Channel.fromPath(params.data_dir))
  //deconvolute(createSignature.out)
  //deconvolute.out.view()
  test.out.view()

}