params.list = "path/to/bwa_index"

params.refdir = '/path/to/refs'
params.dataDir = params.refdir 
params.singleEnd = true
params.threads = 4



Channel.fromPath( file(params.list) )
	.splitCsv(header:false, sep :',')
	.map {row ->
		def ref = row[0]
		return [file (ref)]

	}
	.set {input_bwa_index}

process bwa_index {
	publishDir params.dataDir, mode: 'copy' 

	input:
	file(ref) from input_bwa_index
	

	output:
	file ("${ref}*") into bwa_index_output



	shell:
	'''
	bwa index  !{ref}
	
	'''

}






