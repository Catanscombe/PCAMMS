
params.reads = '*.fastq'
params.pear = '/home/ubuntu/metagenomics/tools/pear-0.9.10-bin-64/pear-0.9.10-bin-64'
params.outdir = 'path/to/outdir'
params.read_path = params.outdir + params.reads
params.dataDir = params.outdir + 'fastqc_output'
params.singleEnd = true

Channel
	.fromFilePairs (params.read_path, size: params.singleEnd ? 1 : 2)
	.set {do_fastqc}

	
	

process do_fastqc  {
	publishDir  params.dataDir ,  mode: 'copy', 
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	
	input:
	set val(name), file(reads) from do_fastqc
	
	output:
	file "*_fastqc*" into fastqc_output

	
	
	shell:
	'''
	fastqc !{name}.fastq
	'''
	
	
	
}

