params.reads = '*_R{1,2}.fastq'
params.outdir = '/home/ubuntu/metagenomics/projects/may22'
params.read_path = params.outdir  + params.reads
params.dataDir = params.outdir + '/neg_map'
params.singleEnd = true
params.threads = 4
//params.contamination = '*_contamination.fasta'
params.contamination_path = 'path/to/contam'

Channel
	.fromFilePairs (params.read_path, flat: true)
	.ifEmpty{exit 1 , 'found no files'}
	.set{reads}


process map_neg {
	//publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file(forward), file(reverse) from reads
	output:
	set dataset_id , file ("${dataset_id}_map_neg.sam") into map_neg_output


	shell:
	'''
	bwa mem -t !{params.threads} !{params.contamination_path} !{forward}  !{reverse} > !{dataset_id}_map_neg.sam

	'''

}


process extract_unmapped {
	publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file("${dataset_id}_map_neg.sam") from map_neg_output

	output:
	file ("${dataset_id}_R1.fastq") 
	file ("${dataset_id}_R2.fastq")

	shell:
	'''
	samtools view -h -b -S !{dataset_id}_map_neg.sam > !{dataset_id}_map_neg.bam
	samtools view -b -f 4 !{dataset_id}_map_neg.bam > !{dataset_id}_unmapped.bam
	bamToFastq -i  !{dataset_id}_unmapped.bam -fq !{dataset_id}_R1.fastq -fq2 !{dataset_id}_R2.fastq 
	'''

}

