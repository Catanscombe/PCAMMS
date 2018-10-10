params.input = "path/to/reference/mapped"
params.output = "path/to/mapped"
params.list = "path/to/list"


Channel
	.fromPath(file(params.list))
	.splitCsv(header:false)
	.map { row ->
		def sample = row[0]
		def organism  = row[2]
		def name = sample + organism
		def input_bam = name + '.bam'

		return[sample, organism, name, file(input_bam)]

	}
	.set{input_extract_mapped}


process extract_mapped {
	publishDir params.output, mode: 'copy'

	input:
	set sample, organism, name, file(input_bam) from input_extract_mapped

	output:
	file("${name}_R1.fastq")
	file("${name}_R2.fastq")


	shell:
	'''
	samtools view -b -F 4 !{params.input}!{input_bam} > !{name}_mapped.bam
	bamToFastq -i !{name}_mapped.bam -fq !{name}_R1.fastq -fq2 !{name}_R2.fastq
	'''



}