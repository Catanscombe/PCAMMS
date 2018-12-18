# PCAMMS
Pipeline for Classifying and Assembling Microbes from Metagenomic Sequencing 

This pipeline is still in development

Description:

Brief:
This pipeline is for initial analysis of metagenomic samples sequenced on the Miseq platforms. It was originally designed for detection of pathogens in sterile site infections. The pipeline consists if two scripts, the first (metagenomic_qc.py)  prepares the reads for analysis (including merging reads, removing host and performing basic QC). The second (metagenomic.py) performs classification on the reads using CLARK in light mode, then downloads the representative reference file for each species taxon ID estimated to be significant, and reference assemblies the unpaired reads against his reference. 
 
In-depth:
*Metagenomic_qc.py*

Will copy the raw files to output directory and rename based either on the sample number and run name provided (default), a given sample list.  Option available to not rename the files, raw files remain unchanged. 

Remove host reads by mapping against the host reference, default is human. These fastq files replace the raw files in the output directory. Host can be changed see below. Option to not map against host.  

Merge paired end reads

Perform basic qc using Fastqc on merged reads, then produce a combined report using mulitqc. 

Classify negative control(s) (optional)


*metagenomic.py*

Optional: Map against the negative control(s) and/or a local contamination library to remove any contaminating reads

Classify the samples using CLARK in light mode

Use the average read length and number of reads classified to assess the possibility of a true hit. (av_readlength x number of classified hits)/genome size if this is greater than 0.001 (0.1%) then the hit is considered significant
This balances not missing hits, but also not performing assembly on very low quality hits. 


Installation:

<<<<<<< HEAD
Please note building the database requires 90GB of storage. Once build is complete raw files are removed and the total required space is 2.2GB


=======
>>>>>>> ca102a577722815324bc12b79ca057f0ff1214ba
git clone https://github.com/Catanscombe/PCAMMS.git

cd PCAMMS

nohup ./setup



this will install the following prerequisites:
*java version 8 

*fastq-join (https://github.com/brwnj/fastq-join)

*CLARK CLAssifier based on Reduced K-mers (http://clark.cs.ucr.edu) 
	*including building database
	
*Krona tools (https://github.com/marbl/Krona/wiki/KronaTools)

*Nextflow (https://www.nextflow.io)

*Fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

*Multiqc (https://multiqc.info)

*seqtk (https://github.com/lh3/seqtk)

*bwa (http://bio-bwa.sourceforge.net)

Usage:
usage: metagenomic_qc.py [-h] [-nl NEG_LIST] [-r {yes,no}] [-id SAMPLE_ID]
                         [-s SAMPLE_LIST] [-n1 NEG_SAMPLE_R1]
                         [-n2 NEG_SAMPLE_R2] [-t THREADS] [-rh {yes,no}]
                         input_dir output_dir run


Basic

*metagenomic_qc.py*

python Path/to/PCAMMS/metagenomic_qc.py directory/containing/raw-illumina/fastqs output/directory runname

outputs:
directory will be made in output/directory/qc_results 
contains:
multiqc_report.html: a webpage containing basic qc data performed on the merged reads
runname_host_map_info.csv: contains details of the number and % of reads which mapped to the host genome
runname_read_info.csv: contains raw read number, merged read number, % reads merged

*metagenomic.py*

python Path/to/PCAMMS/metagenomic.py output/directory runname

outputs:

directory will be made in output/directory/results 
contains:
multiqc_report.html: a webpage containing basic qc data performed on the merged reads

