#!/usr/bin/env nextflow

//create params
//rather Channel.from (multiple values but only use once) or Channel.value (one value, use multiple times)?
//params: single cell dataset, rnaseq dataset, method??
params.single_cell_list = ["Maynard", "test"]
params.rnaseq_list = ["hoek", "finotello", "whatever"]
params.method_list =['bisque', 'cibersortx', 'scaden'] // many more
params.dataFolder = "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/"
params.saveFolder = "/nfs/data/omnideconv_benchmarking/"


process createSignature {

	input:
	val sc from Channel.fromList(params.single_cell_list)
	each rnaseq from Channel.fromList(params.rnaseq_list) // can i use each mode here as well? since mode is already used.. is this necessary?
	each mode from params.method_list
	path data from params.dataFolder
	path saveFolder from params.saveFolder

	output:
	stdout into signatures
	//file *.RData into signatures

	"""
	echo '$sc $data $mode $rnaseq'
	Rscript /nfs/proj/omnideconv_benchmarking/benchmark/pipeline/computeSignaturesNF.R --path '$data' --single_cell '$sc' --rna_seq '$rnaseq' --method '$mode' --saveFolder '$saveFolder'
	# skip few methods without signature here and just return NULL?
	""" 
}

signatures.view{ it }

