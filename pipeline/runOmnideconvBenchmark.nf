#!/usr/bin/env nextflow

//create params
//rather Channel.from (multiple values but only use once) or Channel.value (one value, use multiple times)?
//params: single cell dataset, rnaseq dataset, method??
params.single_cell_list = ["Maynard", ""]
params.rnaseq_list = ["hoek", "finotello"]
params.method_list =['bisque', 'cibersortx', 'scaden'] // many more
params.data = "file path"

process createSignature {
	input:
	val sc from Channel.fromList(params.single_cell_list)
	val rnaseq from Channel.fromList(params.rnaseq_list) // can i use each mode here as well? since mode is already used.. is this necessary?
	each mode from params.method_list
	path sc_data_path
	output:
	file *.RData into signatures

	"""
	// call Rscript for each possible combination?
	// skip few methods without signature here and just return NULL? (if() \n """""" \n else if() \n """""")
	""" 
