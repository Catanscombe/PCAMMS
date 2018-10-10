params.reads = '*_R{1,2}.fastq'
params.input = '/path/to/mapped/'

params.read_pairs = params.input + params.reads

Channel
	.fromFilePairs(params.read_pairs, flat: true)
	.ifEmpty {exit 1, 'Found no input reads'}
	.set {input_denovo}

process denovo {
	tag {pair_id}
	publishDir params.input, mode: 'copy'
	input:
	set pair_id, file(forward), file(reverse) from input_denovo

	output:
	file "${pair_id}" into pear_output

	shell:
	'''
	spades.py -1 !{params.input}!{forward} -2  !{params.input}!{reverse} -o !{pair_id}
	'''

}

