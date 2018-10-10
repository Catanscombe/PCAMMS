params.reads = '*.fastq'
params.outdir = 'path/to/results'
params.read_path = params.outdir  + params.reads
params.dataDir = params.outdir 
params.singleEnd = true
params.threads = 4

params.current_path ='path/to/current'
params.clark_path = 'path/to/clark'

Channel
	.fromFilePairs (params.read_path, size: params.singleEnd ? 1 : 2)
	.ifEmpty{exit 1 , 'found no files'}
	.set{reads}


process classify {
	publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file(reads) from reads
	output:
	set dataset_id , file ("${dataset_id}_clark.csv") into classify_output
	


	shell:
	'''
	
	./classify_metagenome.sh  --light -O !{params.outdir}!{reads} -R !{params.outdir}!{dataset_id}_clark
	
	'''

}

//process abundance {
//	publishDir params.dataDir, mode: 'copy' 
//
//	input:
//	set dataset_id, file("${dataset_id}_clark.csv") from classify_output
//	output:
//	set dataset_id , file ("${dataset_id}_clark_abundance.csv") into abundance_output
//

//	shell:
//	'''
//	cd !{params.clark_path}
//	./estimate_abundance.sh --krona -F !{params.outdir}!{dataset_id}_clark.csv -D DIR_DB > !{params.outdir}!{dataset_id}_clark_abundance.csv
//	cd !{params.current_path}
//	'''

//}


//process krona {
//	publishDir params.dataDir, mode: 'copy' 

//	input:
//	set dataset_id, file("${dataset_id}_clark_abundance.csv") from abundance_output

//	output:
//	set dataset_id , file ("${dataset_id}_clark.html") into krona_output

//	shell:
//	'''
//	ktImportTaxonomy -o !{dataset_id}_clark.html -m 3 !{dataset_id}_clark_abundance.csv
//	'''

//}

