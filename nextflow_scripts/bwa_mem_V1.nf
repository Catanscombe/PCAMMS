params.list = "/path/to/sample_ref.csv"
params.path = "/path/to/paired/"
params.refdir = '/path/to/refs/'
params.dataDir = '/path/to/mapping/'
params.threads = 4


Channel.fromPath( file(params.list) )
		.splitCsv(header:false, sep :',')
		.map {row ->
			def sample = row[0]
			def sample_path_R1 =  sample + '_R1' + '.fastq'
			def sample_path_R2 =  sample + '_R2' + '.fastq'
			def ref = row[6]
			def name = row[1]
			
			return [sample ,file(sample_path_R1), file (sample_path_R2) ,file (ref), name]
			
			
			
		}
		.set {input_bwa_mem}


	
process bwa_mem {
	
	publishDir params.dataDir, mode: 'copy'
	input:
	set sample ,file(sample_path_R1), file(sample_path_R2) , file(ref), name from input_bwa_mem



	output:
	file("${sample}${name}.bam")
	file("${sample}${name}.txt")
	file("${sample}${name}_mapped.txt") 

	shell:
	'''
	bwa mem -t !{params.threads} !{params.refdir}!{ref}  !{params.path}!{sample_path_R1}  !{params.path}!{sample_path_R2} > !{sample}!{name}.sam
	samtools view -h -b -S !{sample}!{name}.sam > !{sample}!{name}_.bam
	samtools sort !{sample}!{name}_.bam > !{sample}!{name}.bam
	samtools depth !{sample}!{name}.bam | wc -l > !{sample}!{name}.txt
	samtools view -F 0x904 -c !{sample}!{name}.bam > !{sample}!{name}_mapped.txt
	'''

}


