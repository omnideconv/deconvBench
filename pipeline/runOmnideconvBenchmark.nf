#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process createSignature {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each sc
	each rnaseq
	each mode 
	path rnaseq_data
	path sc_data

	output:
	tuple val(sc), val(rnaseq), val(mode), file('*.rds'), emit: signature

	"""
	echo '$sc $sc_data $rnaseq_data $mode $rnaseq'
	#./computeSignaturesNF.R --sc_path '$sc_data' --single_cell '$sc' --rnaseq_path '$rnaseq_data' --rna_seq '$rnaseq' --method '$mode'
	Rscript /nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/computeSignaturesNF.R --sc_path '$sc_data' --single_cell '$sc' --rnaseq_path '$rnaseq_data' --rna_seq '$rnaseq' --method '$mode'
	""" 
}

process deconvolute { 

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	tuple val(sc), val(rnaseq), val(mode), val(signature) 
	path rnaseq_data
	path sc_data
	
	output:
	path '*.rds', emit: deconvolution

	"""
	echo '$sc $rnaseq $signature $sc_data $rnaseq_data'
	Rscript /nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/runDeconvolutionNF.R --sc_path '$sc_data' --single_cell '$sc' --rnaseq_path '$rnaseq_data' --rna_seq '$rnaseq' --method '$mode' --signature '$signature'
	""" 
}


workflow{
  single_cell = Channel.fromList(params.single_cell_list)
  rna_seq = Channel.fromList(params.rnaseq_list)
  methods = params.method_list
  rnaseq_data = Channel.fromPath(params.data_dir_rnaseq).collect()
  sc_data = Channel.fromPath(params.data_dir_sc).collect()
  createSignature(single_cell, rna_seq, methods, rnaseq_data, sc_data)
  deconvolute(createSignature.out.signature, rnaseq_data, sc_data)
  deconvolute.out.deconvolution.view()
  
}