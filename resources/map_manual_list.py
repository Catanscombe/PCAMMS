import os                                                                                                                           
import csv
import argparse
import fileinput
 
  
def main():
 	parser = argparse.ArgumentParser()
    parser.add_argument('list', help = 'list containing: sample_ID, ref_ID, ref_file_name.fasta')
    parser.add_argument('reference directory' , help = 'path to directory containing reference files')
    parser.add_argument('results_dir' , help = 'path to results directory from metagenomic.py')
    
    
    args = parser.parse_args()
    ref_mapping(args)
    stats(args)
    make_csv(args, ID)
  	extract_map(ID)
  	denovo(ID)

def ref_mapping(args):
    os.system('bwa index %s ' % (args.reference))
    os.system('bwa mem %s %s %s > %s_%s.sam' % (args.reference, args.forward, args.reverse, args.sample, args.identifier))
  
def stats(args):
    ID = ('%s_%s' % (args.sample, args.identifier))
    os.system('samtools view -h -b -S %s.sam > %s.bam' % (ID, ID))
    os.system('samtools sort %s.bam > %s_s.bam' % (ID, ID))
    os.system('samtools depth %s_s.bam | wc -l > %s.txt' % (ID, ID))
    os.system('samtools view -F 0x904 -c %s_s.bam > %s.mapped.txt' % (ID, ID))
 	
 	os.system('rm %s.sam %s.bam ' % (ID, ID))

 
def make_csv(args, ID):
    ID = ('%s_%s' % (args.sample, args.identifier))
    os.system ("grep -v '>' %s | wc | awk '{print $3-$1}' > %s_len.txt" % (args.reference, args.reference))
 
    with open ('%s_ref_map.csv' % (ID), 'a') as f:

	    for line in open ('%s.txt' % (ID)):
            line =line.strip()
            cov = float(line)
  
        for line in open ('%s.mapped.txt' % (ID)):
            line = line.strip()
            reads = float(line)
  
        for line in open ('%s_len.txt' % (args.reference)):
            line = line.strip()
            ref_len = float(line)
  		if cov > 0:
        	cov_pec = cov/ref_len
        else:
        	cov_pec = 0 
        writer = csv.writer(f, delimiter = ',')
        writer.writerow([ID, reads, cov, ref_len, cov_pec])
 
    for line in fileinput.input(files = ['%s_ref_map.csv' % (ID)], inplace =True):
        if fileinput.isfirstline():
            print 'ID, reads mapped, cov , ref_len, % cov'
        print line,

def extract_map(ID):
	os.system ('samtools view -b -F 4 %s_s.bam  > %s_map.bam ' % (ID, ID))
	os.system ('bamToFastq -i %s_map.bam -fq %s.fq  ' % (ID, ID))
  
def denovo(ID):
	 os.system('spades.py -s %s.fq  -o %s' % (ID, ID))
	 try:
	 	os.system ('cp %s/contigs.fasta %s_spades.fasta' % (ID, ID))
	 else:
	 	pass  
main()

