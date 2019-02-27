import os                                                                                                                           
import csv
import argparse
import fileinput
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('list', help = 'list of samples to extract, same format as reference_assembled.csv')
    parser.add_argument('results_dir' , help = 'path to results directory from metagenomics.py')
    args = parser.parse_args()
    ID = extract_map(args)
    denovo(ID)
 
def extract_map(args):


	with open ('%s' % (args.list), 'r') as f:
		next(f)
		for line in f:
			line = line.strip().split(',')
			sample = line[0]
			identifier = line[2]
			taxon_ID = [8]

			ID = ('%s%s' % (sample, identifier))

    
    		os.system ('samtools view -b -F 4 %s/reference_mapping/%s.bam  > %s_map.bam ' % (args.results_dir, ID, ID))
    		os.system ('bamToFastq -i %s_map.bam -fq %s.fq ' % (ID, ID))


			os.system('spades.py -s %s.fq -o %s' % (ID, ID))
    		os.system ('cp %s/contigs.fasta %s_spades.fasta' % (ID, ID))

main()