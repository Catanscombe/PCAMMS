import os                                                                                                                           
import csv
import argparse
import fileinput
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', help = 'sample_ID')
    parser.add_argument('identifier' , help = 'string to identify this denovo assembly')
    parser.add_argument('bam' , help = 'path to bam from reference mapping (found in run/reference_mapping directory)')
    args = parser.parse_args()
    ID = extract_map(args)
    denovo(ID)
 
def extract_map(args):
    ID = ('%s_%s' % (args.sample, args.identifier))
    os.system ('samtools view -b -F 4 %s  > %s_map.bam ' % (args.bam, ID))
    os.system ('bamToFastq -i %s_map.bam -fq %s.fq ' % (ID, ID))
 
    return ID
 
def denovo(ID):
    os.system('spades.py -s %s.fq -o %s' % (ID, ID))
    os.system ('cp %s/contigs.fasta %s_spades.fasta' % (ID, ID))

main()