import os                                                                                                                           
import csv
import argparse
import fileinput


def main():
 	parser = argparse.ArgumentParser()
    parser.add_argument('sample', help = 'sample_ID')
    parser.add_argument('fastq' , help = 'path to rinterweaved fastq')
    parser.add_argument('clark', help = 'path to clark csv file (found in sample/sample_clark.csv)')
    parser.add_argument('identifier' , help = 'string to identify this refernce assembly')
    parser.add_argument('taxon_ID', help =  'taxon ID')

    args = parser.parse_args()

def extract(args):

	ID = ('%s_%s' % (args.sample, args.identifier))

	os.system ('grep %s %s | awk -F , "{print $1}" > %s_clark.txt' %  (args.taxon_ID, args.clark, ID) )

	os.system ('seqtk subseq %s %s_clark.txt > %s.fq' % (args.fastq, ID, ID))

main()
