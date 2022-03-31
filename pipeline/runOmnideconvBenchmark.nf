#!/usr/bin/env nextflow
nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

def check_samplesheet(input_file, base_dir) {
    def samples = []
    for(line in parseCsv(new FileReader(input_file))) {
        expression_matrix = file("${base_dir}/${line['matrix_path']}")
        meta = line.toMap()
        meta.remove('matrix_path')
        samples << [meta, expression_matrix]
    }
    return samples
}

def check_samplesheet_summarized(input_file, base_dir) {
    def samples = [:]
    for(line in parseCsv(new FileReader(input_file))) {
        expression_matrix = file("${base_dir}/${line['matrix_path']}")
        set = line['dataset']
        type = line['type']
        println "$set $type"
        if(samples.containsKey(set)){
          samples[set].add(set + "_" + type + ":" + expression_matrix)
        } else {
          samples[set] = []
          samples[set].add(set + "_" + type + ":" + expression_matrix)
        }
    }
    return samples.values()
}

process setup {

	publishDir params.results_dir, mode: params.publishDirMode

	output:
	path('*.tsv'), emit: runtime

	"""
	printf "method\tsc_dataset\trnaseq_dataset\tstep\truntime_user\truntime_system\truntime_elapsed\n" > runtimes.tsv
	""" 
}


process createSignature {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each rnaseq
	each mode 
	path rnaseq_dir
	tuple val(sc_meta), path(sc_matrix)
	path sc_dir
	path mapping_sheet

	output:
	tuple val(sc_meta), path(sc_matrix), val(rnaseq), val(mode), path('*.rds'), emit: signature
	stdout emit: runtime

	shell:
	'''
	#echo '$sc $sc_dir $rnaseq_dir $mode $rnaseq'
	#computeSignaturesNF.R '$sc_dir' '$sc' '$rnaseq_dir' '$rnaseq' '$mode' '$mapping_sheet'
	TIME="$(computeSignaturesNF.R '!{sc_matrix}' '!{sc_meta}' '!{sc_dir}' '!{rnaseq_dir}' '!{rnaseq}' '!{mode}' '!{mapping_sheet}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc_meta}'\t'!{rnaseq}'\tbuild_model\t'$FINAL
	''' 
}

process deconvolute { 

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	tuple val(sc_meta), path(sc_matrix), val(rnaseq), val(mode), val(signature)
	path rnaseq_dir
	path sc_dir
	path mapping_sheet

	output:
	path('*.rds'), emit: deconvolution
	stdout emit: runtime

	shell:
	'''
	#echo '$sc $rnaseq $signature $sc_dir $rnaseq_dir'
	#runDeconvolutionNF.R '$sc_dir' '$sc' '$rnaseq_dir' '$rnaseq' '$mode' '$signature' '$mapping_sheet'
	TIME="$(runDeconvolutionNF.R '!{sc_matrix}' '!{sc_meta}' '!{sc_dir}' '!{rnaseq_dir}' '!{rnaseq}' '!{mode}' '!{signature}' '!{mapping_sheet}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc_meta}'\t'!{rnaseq}'\tdeconvolute\t'$FINAL
	''' 
}

process benchmarkingPlots {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	val deconvolution
	path facs
	path mapping_sheet

	
	output:
	path '*.jpeg', emit: plotjpeg

	"""
	plotGtruthHoekNF.R '$deconvolution' '$facs' '$mapping_sheet'
	""" 
}

process simulateBulk {

	errorStrategy 'ignore'
	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each sc_meta
	each scenario 
	path sc_dir
	path mapping_sheet
	val celltypes

	output:
	path('*.rds'), emit: bulk
	//tuple path('*.rds'), val(sc), emit: bulk_tuple

	//path('*.pdf'), emit: correlation

	"""
	echo '$sc_meta $sc_dir $scenario $celltypes'
	simulateBulkNF.R '$sc_dir' '$sc_meta' '$mapping_sheet' '$scenario' 'blood' 'Homo sapiens' '$celltypes'
	""" 
}

process deconvoluteSimulatedBulk {
	
	//errorStrategy 'ignore'
	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each method
	each bulk
	path sc_dir
	val meta
	path mapping_sheet

	output:
	path('*.rds'), emit: deconvolution
	path('*.pdf'), emit: correlation

	"""
	#echo '$sc_dir $method $bulk'
	echo $meta
	deconvoluteSimBulkNF.R '$sc_dir' '$meta' '$bulk' '$method' '$mapping_sheet'
	""" 
}

process plotSimulation {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	val deconv
	val bulk
	path mapping_sheet
	val type

	output:
	path('*.pdf'), emit: entry
	path('*.jpeg'), emit: entry_jpeg

	"""
	plotSimulationCorrelationNF.R '$deconv' '$bulk' '$mapping_sheet' '$type'
	""" 
}

process plotSimulationSpillover {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	val deconv

	output:
	path('*.pdf'), emit: chords

	"""
	plotSpilloverNF.R '$deconv'
	""" 
}


workflow{
  //single_cell = Channel.fromList(params.single_cell_list)
  rnaseq_sets = Channel.fromList(params.bulk_list)
  methods = params.method_list
  scenarios = params.simulation_scenarios
  rnaseq_dir = Channel.fromPath(params.data_dir_bulk).collect()
  sc_dir = Channel.fromPath(params.data_dir_sc).collect()
  remapping = Channel.fromPath(params.mapping_sheet).collect()
  ch_input_files = check_samplesheet("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/sc_files_table_tpmMaynard.csv", params.data_dir_sc)
  Channel.from(ch_input_files).view()
  createSignature(rnaseq_sets, methods, rnaseq_dir, Channel.from(ch_input_files), sc_dir, remapping)
  signatureRuntimes = createSignature.out.runtime
  deconvolute(createSignature.out.signature, rnaseq_dir, sc_dir, remapping)
  deconv = deconvolute.out.deconvolution
  deconvolute.out.deconvolution.view()
  signatureRuntimes.join(deconvolute.out.runtime).collectFile(name: "runtimesCPM.tsv", newLine: true)
  //hoek_samples = deconv.filter{ ds -> ds =~/_hoek/ }                      das
  //hoek_samples_list = hoek_samples.toList()                               hier
  //benchmarkingPlots(hoek_samples_list, rnaseq_dir, rnaseq_sets, remapping)   nicht
  //benchmarkingPlots(deconv.toList(), rnaseq_dir, remapping)
  //benchmarkingPlots.out.plotjpeg.view()
  celltypes = Channel.value("B cell,T cell regulatory (Tregs),T cell CD8+,T cell CD4+,Monocyte non-conventional,Macrophage,Monocyte conventional,NK cell,Myeloid dendritic cell,Plasmacytoid dendritic cell")
  ch_input_files_summarized = check_samplesheet_summarized("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/sc_files_table.csv", params.data_dir_sc)
  //simulateBulk(Channel.from(ch_input_files_summarized), scenarios, sc_dir, remapping, celltypes)
  //sims = simulateBulk.out.bulk.collect()
  //deconvoluteSimulatedBulk(methods, sims, sc_dir, Channel.from(ch_input_files_summarized).collect(), remapping)
  //deconvSim = deconvoluteSimulatedBulk.out.deconvolution
  //true_fractions = deconvSim.filter{ str -> str =~/_truefractions/ }.collect()
  //true_fractions.view()
  //sims.view()
  //plotSimulation(true_fractions, sims, remapping, "truefractions")
  //plotSimulation.out.entry.view()
  //uniform = deconvSim.filter{ str -> str =~/_uniform/ }.collect()view()
  //plotSimulation(uniform, sims, remapping)
  //plotSimulation.out.entry.view()
  //spillover = deconvSim.filter{ str -> str =~/_spillover/ }.collect()
  //plotSimulationSpillover(spillover)
}

