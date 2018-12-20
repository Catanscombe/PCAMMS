params.reads = '*.fna'
params.refdir = '/path/to/refs'
params.read_path = params.refdir  + '/' + params.reads
params.dataDir = params.refdir 
params.singleEnd = true
params.threads = 4



Channel
	.fromFilePairs (params.read_path, size: params.singleEnd ? 1 : 2)
	.ifEmpty{exit 1 , 'found no files'}
	.set{reads}


process bwa_index {
	publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file(reads) from reads
	output:
	set dataset_id , file ("${dataset_id}.fna*") into bwa_index_output



	shell:
	'''
	if [ ! -f !{dataset_id}.fna.amb ]; 
	then 
		bwa index  !{dataset_id}.fna
	fi
	'''

}






