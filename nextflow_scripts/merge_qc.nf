params.reads = '*_R{1,2}.fastq'
params.outdir = '/home/ubuntu/metagenomics/projects/testing/'
params.read_pairs = params.outdir + 'paired/' + params.reads

Channel
	.fromFilePairs(params.read_pairs, flat: true)
	.ifEmpty {exit 1, 'Found no input reads'}
	.set {input_pear}

process pear {
	tag {pair_id}

	input:
	set pair_id, file(forward), file(reverse) from input_pear

	output:
	file "${pair_id}.fastq" into pear_output

	shell:
	'''
	fastq-join  !{forward}  !{reverse} -o !{pair_id}
	mv !{pair_id}join !{pair_id}.fastq
	
	'''

}
pear_output.subscribe {it.copyTo(params.outdir)}

