import os                                                                                                                           
import csv
import argparse
import fileinput


def main():
 	parser = argparse.ArgumentParser()
    parser.add_argument('list', help = 'list of samples to extract, same format as reference_assembled.csv')
    parser.add_argument('results_dir' , help = 'path to results directory from metagenomics.py')


    args = parser.parse_args()
    extract(args)
    

def extract(args):

	with open ('%s' % (args.list), 'r') as f:
		next(f)
		for line in f:
			line = line.strip().split(',')
			sample = line[0]
			identifier = line[2]
			taxon_ID = line[8]

			ID = ('%s%s' % (sample, identifier))

			os.system ("grep %s %s/%s/%s_clark.csv | awk -F ',' '{print $1}'' > %s_clark.txt" %  (taxon_ID, args.results_dir, sample, sample, ID) )

			os.system ('seqtk subseq %s/%s.fastq %s_clark.txt > %s.fq' % (args.results_dir, sample , ID, ID))

			os.system('spades.py -s %s.fq  -o %s' % (ID, ID))
    		os.system ('cp %s/contigs.fasta %s_spades.fasta' % (ID, ID))
	    	
 

    
main()



