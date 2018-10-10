params.reads = '*.fastq'
params.outdir = '/home/ubuntu/metagenomics/projects/may22'
params.read_path = params.outdir  + params.reads
params.dataDir = params.outdir + '/neg_map'
params.singleEnd = true
params.threads = 4
//params.contamination = '*_contamination.fasta'
params.contamination_path = 'path/to/contam'

Channel
	.fromFilePairs (params.read_path, size: params.singleEnd ? 1 : 2)
	.ifEmpty{exit 1 , 'found no files'}
	.set{reads}


process map_neg {
	//publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file(reads) from reads
	output:
	set dataset_id , file ("${dataset_id}_map_neg.sam") into map_neg_output


	shell:
	'''
	bwa mem -t !{params.threads} !{params.contamination_path} !{dataset_id}.fastq > !{dataset_id}_map_neg.sam

	'''

}


process extract_unmapped {
	publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file("${dataset_id}_map_neg.sam") from map_neg_output

	output:
	set file ("${dataset_id}.fastq") into bam_output

	shell:
	'''
	samtools view -h -b -S !{dataset_id}_map_neg.sam > !{dataset_id}_map_neg.bam
	samtools view -b -f 4 !{dataset_id}_map_neg.bam > !{dataset_id}_unmapped.bam
	bamToFastq -i  !{dataset_id}_unmapped.bam -fq !{dataset_id}.fastq  
	'''

}

