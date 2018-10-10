params.reads = '*_R{1,2}.fastq'
params.pear = '/home/ubuntu/metagenomics/tools/pear-0.9.10-bin-64/pear-0.9.10-bin-64'
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
	!{params.pear} -f !{forward} -r !{reverse} -o !{pair_id}
	mv !{pair_id}.assembled.fastq !{pair_id}.fastq
	
	'''

}
pear_output.subscribe {it.copyTo(params.outdir)}

