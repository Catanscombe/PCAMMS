params.reads = '*_R{1,2}.fastq'
params.outdir = '/path/to/reads/'
params.read_path = params.outdir  + params.reads
params.dataDir = params.outdir + '/map_host'
params.singleEnd = true
params.threads = 4

params.host_path = 'path/to/host'

Channel
	.fromFilePairs (params.read_path, flat:true)
	.ifEmpty{exit 1 , 'found no files'}
	.set{reads}


process map_human {
	//publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file(forward), file(reverse) from reads
	output:
	set dataset_id , file ("${dataset_id}_map_host.sam") into map_neg_output


	shell:
	'''
	bwa mem -t !{params.threads} !{params.host_path} !{forward} !{reverse} > !{dataset_id}_map_host.sam

	'''

}


process extract_unmapped {
	publishDir params.dataDir, mode: 'copy' 

	input:
	set dataset_id, file("${dataset_id}_map_host.sam") from map_neg_output

	output:
	file ("${dataset_id}_R1.fastq") into fastq_output
	file ("${dataset_id}_R2.fastq") into fastq_output2_output

	shell:
	'''
	samtools view -h -b -S !{dataset_id}_map_host.sam > !{dataset_id}_map_host.bam
	samtools view -b -f 4 !{dataset_id}_map_host.bam > !{dataset_id}_host_unmapped.bam
	bamToFastq -i  !{dataset_id}_host_unmapped.bam -fq !{dataset_id}_R1.fastq -fq2 !{dataset_id}_R2.fastq  
	'''

}

