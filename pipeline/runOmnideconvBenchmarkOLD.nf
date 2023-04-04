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

def check_samplesheet_individual(input_file, base_dir) {
    def sample = []
    for(line in parseCsv(new FileReader(input_file))) {
		//if(params.single_cell_norm.contains(line['type'])){
		if(params.single_cell_list.contains(line['dataset']) && params.single_cell_norm.contains(line['type'])){	
			//$line.type == params.single_cell_norm
			expression_matrix = file("${base_dir}/${line['matrix_path']}")
        	meta = line.toMap()
        	meta.remove('matrix_path')
        	sample << [meta, expression_matrix]	
		}

    }
	//println sample
    return sample
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
	printf "method\tsc_dataset\trnaseq_dataset\trnaseq_norm\tstep\truntime_user\truntime_system\truntime_elapsed\n" > runtimes.tsv
	""" 
}

process preProcessSingleCell {

	publishDir params.preProcess_dir, mode: params.publishDirMode

	input: 
	tuple val(sc_meta), path(sc_matrix)
	path sc_dir
	each ct_percentage
	each replicates
	val preProcess_dir

	output:
	tuple val(sc_meta), path(sc_matrix), val(rnaseq), val(mode), val(rnaseq_norm), val(results_dir_gen), emit: dataset

	shell:
	'''
	/vol/omnideconv/benchmark/pipeline/bin/preProcessSingleCell.R '!{sc_matrix}' '!{sc_meta}' '!{sc_dir}' '!{ct_percentage}' '!{replicates}' '!{preProcess_dir}' '!{seed}'
	''' 
}


process createSignature {
	label "deconvolution_related"
	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each rnaseq
	each mode 
	path rnaseq_dir
	each rnaseq_norm
	tuple val(sc_meta), path(sc_matrix)
	path sc_dir
	val results_dir_gen


	output:
	tuple val(sc_meta), path(sc_matrix), val(rnaseq), val(mode), val(rnaseq_norm), val(results_dir_gen), emit: signature
	stdout emit: runtime

	shell:
	'''
	TIME="$(/vol/omnideconv/benchmark/pipeline/bin/computeSignaturesNF.R '!{sc_matrix}' '!{sc_meta}' '!{sc_dir}' '!{rnaseq_dir}' '!{rnaseq}' '!{rnaseq_norm}' '!{mode}' '!{results_dir_gen}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc_meta}'\t'!{rnaseq}'\tbuild_model\t'$FINAL
	''' 
}

process deconvolute { 
	label "deconvolution_related"
	publishDir params.results_dir, mode: params.publishDirMode

	input:
	tuple val(sc_meta), path(sc_matrix), val(rnaseq), val(mode), val(rnaseq_norm), val(results_dir_gen)
	path rnaseq_dir
	path sc_dir
	

	output:
	tuple val(sc_meta), val(rnaseq), val(rnaseq_norm), val(mode), val(results_dir_gen), emit: deconvolution
	stdout emit: runtime

	shell:
	'''
	TIME="$(/vol/omnideconv/benchmark/pipeline/bin/runDeconvolutionNF.R '!{sc_matrix}' '!{sc_meta}' '!{sc_dir}' '!{rnaseq_dir}' '!{rnaseq}' '!{rnaseq_norm}' '!{mode}' '!{results_dir_gen}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc_meta}'\t'!{rnaseq}'\tdeconvolute\t'$FINAL
	''' 
}

process computeMetrics {
	label "deconvolution_related"
	publishDir params.results_dir, mode: params.publishDirMode

	input:
	tuple val(sc_meta),  val(rnaseq), val(rnaseq_norm), val(mode), val(results_dir)
	val rnaseq_dir

	output:
	stdout emit: confirm_execution

	shell:
	'''
	echo '!{sc_meta}' '$rnaseq $rnaseq_norm $mode $results_dir'
	/vol/omnideconv/benchmark/pipeline/bin/computeMetricsNF.R '!{sc_meta}' '!{rnaseq}' '!{rnaseq_norm}' '!{mode}' '!{results_dir}' '!{rnaseq_dir}'
	echo 'run completed'
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
	label "deconvolution_related"
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
	label "deconvolution_related"
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
  //setup()
  rnaseq_sets = Channel.fromList(params.bulk_list)
  methods = params.method_list
  rnaseq_norm = params.bulk_norm
  scenarios = params.simulation_scenarios
  rnaseq_dir = Channel.fromPath(params.data_dir_bulk).collect()
  sc_dir = Channel.fromPath(params.data_dir_sc).collect()
  remapping = Channel.fromPath(params.mapping_sheet).collect()
  ch_input_files = check_samplesheet_individual(params.sc_files, params.data_dir_sc)
  Channel.from(ch_input_files).view()
  createSignature(rnaseq_sets, methods, rnaseq_dir, rnaseq_norm, Channel.from(ch_input_files), sc_dir, params.results_dir_general)
  signatureRuntimes = createSignature.out.runtime
  deconvolute(createSignature.out.signature, rnaseq_dir, sc_dir)
  deconv = deconvolute.out.deconvolution
  deconvolute.out.deconvolution.view()
  signatureRuntimes.mix(deconvolute.out.runtime).collectFile(storeDir: params.results_dir, name: file("runtimes.tsv"), newLine: true)
  //computeMetrics(deconvolute.out.deconvolution, rnaseq_dir)
  //benchmarkingPlots(deconv.toList(), rnaseq_dir, remapping)
  //benchmarkingPlots.out.plotjpeg.view()
  //celltypes = Channel.value("B cell,T cell regulatory (Tregs),T cell CD8+,T cell CD4+,Monocyte non-conventional,Macrophage,Monocyte conventional,NK cell,Myeloid dendritic cell,Plasmacytoid dendritic cell")
  //ch_input_files_summarized = check_samplesheet_summarized("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/sc_files_table.csv", params.data_dir_sc)
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