import argparse
import os
import csv
import sys
from collections import defaultdict
import fileinput
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', help ='output directory')
    parser.add_argument( 'run',  help = 'run description')
    parser.add_argument ('-c' , '--contamination' , help = 'path to fasta containing local contamination library')
    parser.add_argument('-nc' , '--neg_control' , nargs='+' ,  help = 'list the merged negative control fastq files from metagenomic.py you wish to be used as contamination libraries')
    parser.add_argument('-tc' , '--taxon_contaminants' , help = 'provide a list of taxon ID you wish to be ignored in the process')
    parser.add_argument ('-t' , '--threads' , default = '4', help = 'number of threads availble, default is 4')
    

    args = parser.parse_args()
    
    directory = os.path.dirname(os.path.abspath(sys.argv[0]))
   
    
    check(args)
    
    neg_library(args, directory)
    results_dir = classify_samples(args, directory)
    classification_info(args, results_dir)
    predict_genome_cov(args, directory)
    auto_assemble(args, directory)
    taxon_abundance(args)
    remove_work_dir ()

def check(args):
# check the input parameters
    counter = 0

    if os.path.isdir(args.output_dir):
        print 'output_dir found'
        if not any(name.endswith('.fastq') for name in os.listdir(args.output_dir)):
            print 'directory has no fastq files, please check path'
            counter = counter +1
    else:
        print 'cannot access output_dir, please check path'
        counter = counter +1

    if args.contamination:
        if os.path.isfile(args.contamination):
            print 'local contamination library given'

            if args.contamination.endswith('.fa' or '.fasta'):
                print 'local contamination file given in fasta format'
            else:
                print 'please supply local contamiation library in fasta format'
                counter = counter +1 
        else:   
            print 'cannot access contamination library'

    if args.neg_control:
        print args.neg_control
        for x in args.neg_control:
            if os.path.isfile(x):
                print ('%s found, will be used as a contamination library' % (x))
                if x.endswith('.fastq'):
                    print ('fastq file entered')
                else:
                    print ('fastq expected')
                    counter = counter +1
            else:
                print ('%s connot be found please check path' % (x))
                counter = counter +1
     
    if args.taxon_contaminants:
        if os.path.isfile(args.taxon_contaminants):
            print ('list of contamination taxon ID provided, taxon listed in %s will be ignored from the analysis' % (args.taxon_contaminants))

    if counter >= 1:
        print ('problems with %s  input file or directory' %(counter))
        sys.exit()



def neg_library(args, directory):
    if args.contamination:
        try:
            os.mkdir('%s/contamination' % (args.output_dir))
        except OSError:
            pass

        os.system ('cp %s %s/contamination' % (args.contamination, args.output_dir))

    if args.neg_control:
        try:
            os.mkdir('%s/contamination' % (args.output_dir))
        except OSError:
            pass

        for x in args.neg_control:
            x_split = x.split('/')
            x_file = x_split[-1]
            x_file_split = x_file.split('.')
            x_id = x_file_split[0]
 
            os.system('seqtk seq -a %s > %s/contamination/%s.fasta' % (x,args.output_dir, x_id))

    if args.contamination or args.neg_control:
        os.system ('cat %s/contamination/*.fasta > %s/contamination/%s_contamination.fasta' % (args.output_dir, args.output_dir, args.run))
        for each_file in os.listdir('%s/contamination' % (args.output_dir)):
            if each_file != ('%s_contamination.fasta' % (args.run)):
                os.system('rm %s/contamination/%s' % (args.output_dir, each_file))

        os.system('bwa index %s/contamination/%s_contamination.fasta' % (args.output_dir, args.run))
        contam = os.path.abspath('%s/contamination/%s_contamination.fasta' % (args.output_dir, args.run))

        try:
            os.mkdir('%s/neg_map' % (args.output_dir))
        except OSError:
            pass
        if args.threads < 4:
            os.system ('nextflow run %s/nextflow_scripts/neg_map.nf --contamination_path "%s" --outdir "%s/" --threads "%s"' % (directory, contam, args.output_dir, args.threads))

        else: 
            os.system ('nextflow run %s/nextflow_scripts/neg_map.nf --contamination_path "%s"  --outdir "%s/" ' % (directory, contam, args.output_dir))

        try:
            os.mkdir('%s/results' % (args.output_dir))
        except OSError:
            pass
        
        with open ('%s/results/%s_contamination_map_info.csv' % (args.output_dir, args.run), 'a') as f:
            writer = csv.writer (f, delimiter = ',')
            for each_file in os.listdir ('%s' % (args.output_dir)):
                if each_file.endswith('.fastq'):
                    each_file_split = each_file.split('.')
                    sample = each_file_split[0]
                    merge_count = 0
                    for line in open ('%s/%s' % (args.output_dir, each_file)).xreadlines( ): merge_count +=1
                    neg_unmap_count = 0
                    for line in open ('%s/neg_map/%s.fastq' % (args.output_dir, sample)).xreadlines(): neg_unmap_count +=1
                    merge_reads = merge_count/4
                    neg_unmap_reads = neg_unmap_count/4
                    neg_map_reads =  merge_reads - neg_unmap_reads
                    pec_neg_map = float(neg_map_reads)/float(merge_reads)
                    writer.writerow([sample, merge_reads, neg_map_reads , pec_neg_map, neg_unmap_reads])
        
        os.system('sort %s/results/%s_contamination_map_info.csv > %s/results/%s_contamination_map_info2.csv' % (args.output_dir, args.run, args.output_dir, args.run))
        os.system ('mv %s/results/%s_contamination_map_info2.csv  %s/results/%s_contamination_map_info.csv' % (args.output_dir, args.run, args.output_dir, args.run))
        for line in fileinput.input(files=['%s/results/%s_contamination_map_info.csv' % (args.output_dir, args.run)] , inplace = True):
            if fileinput.isfirstline():
                print 'sample,input_reads,contamination_map_reads,pec_contamination_map,output_reads'
            print line,
        os.system('cp %s/neg_map/*.fastq %s/' % (args.output_dir, args.output_dir))
        os.system('rm -r  %s/neg_map' % (args.output_dir))

        print 'contamination mapping complete, contamination reads removed from further analysis'

def classify_samples(args, directory):
    Clark_dir = ('%s/CLARKSCV1.2.5.1' % (directory))
    
    for each_file in os.listdir(args.output_dir):
        if each_file.endswith('.fastq'):
            each_file_split = each_file.split('.')
            sample = each_file_split[0]
            try:
                os.mkdir('%s/%s' % ( args.output_dir, sample))
            except OSError:
                pass
            cwd = os.getcwd()
            results_dir = os.path.dirname(os.path.abspath('%s/%s' % (args.output_dir, sample)))
            print results_dir
            os.chdir('%s' % (Clark_dir))
            os.system ( './classify_metagenome.sh  --light -O %s/%s.fastq  -n %s -R %s/%s/%s_clark' % ( results_dir, sample, args.threads , results_dir, sample, sample))
            os.system ('./estimate_abundance.sh --krona -F %s/%s/%s_clark.csv -D DIR_DB' % ( results_dir, sample, sample))
            os.system ('mv results.krn %s/%s/%s_clark_abundance.csv' % (results_dir, sample, sample))
            os.chdir ('%s' % (cwd))
            os.system  ('ktImportTaxonomy -o  %s/%s/%s.html -m 3 %s/%s/%s_clark_abundance.csv' % ( results_dir, sample,  sample, results_dir, sample, sample))

            print ('classification performed on %s' % (sample))

    print ('classification performed on all samples in %s' % (results_dir))

    return results_dir


def classification_info(args, results_dir):
    try:
        os.mkdir('%s/results' % (args.output_dir))
    except OSError:
        pass

    d_raw = defaultdict(str)
    for line in open ('%s/qc_results/%s_host_map_info.csv' % (results_dir, args.run)).readlines():
        try: 
            line = line.strip().split(',')
            key,value = line[0] , line[1]
            d_raw[key] += value
        except IndexError:
            key, value = 'null', 'null'
            d_raw[key] += value

    


    d_merge = defaultdict(str)
    for line in open ('%s/qc_results/%s_read_info.csv' % (results_dir, args.run)).readlines():
        try:
            line = line.strip().split(',')
            print line[0]
            key, value = line[0], line[2]
            d_merge[key] += value
        except IndexError:
            key, value = 'null' , 'null'
            d_merge[key] += value

    print d_merge

    d_map_host = defaultdict(str)
    for line in open ('%s/qc_results/%s_host_map_info.csv' % (results_dir, args.run)).readlines():
        try:
            line = line.strip().split(',')
            print line[0]
            key, value = line[0], line[2]
            d_map_host[key] += value
        except IndexError:
            key, value = 'null' , 'null'
            d_map_host[key] += value

 
    if args.contamination or args.neg_control:

        d_map_contamination = defaultdict(str)
        for line in open ('%s/results/%s_contamination_map_info.csv' % (results_dir, args.run)).readlines():
            try:
                line = line.strip().split(',')
                print line[0]
                key, value = line[0], line[2]
                d_map_contamination[key] += value
            except IndexError:
                key, value = 'null' , 'null'
                d_map_contamination[key] += value

        print d_map_contamination


    for fastq in os.listdir(args.output_dir):
        if fastq.endswith('.fastq'):
 
            fastq_split = fastq.split('.')
            sample = fastq_split[0]

            raw_reads= int(d_raw[sample])
            merge_reads = int(d_merge[sample])            
            host_map = int(d_map_host[sample])
            host_unmap = raw_reads - host_map
            if host_map > 0 and raw_reads > 0: 
                host_pec = float(host_map)/float(raw_reads)
            else: 
                host_pec  = 0

            if merge_reads > 0 and host_unmap > 0:
                merge_pec = float(merge_reads)/float(host_unmap)
            else:
                merge_pec = 0

            if args.contamination or args.neg_control:
                contamination_map = int(d_map_contamination[sample])
                contamination_unmap = host_unmap - contamination_map
                if contamination_map > 0 and host_map > 0:
                    contamination_pec = float(contamination_map)/float(host_unmap)
                else:
                    contamination_pec = 0 


            unclassified = open('%s/%s/%s_clark.csv' % (results_dir,  sample, sample)).read()
            unclassified_count = unclassified.count('NA')                                   
            
            if args.contamination or args.neg_control:
                classification = int(contamination_unmap) - int(unclassified_count)
                if unclassified_count > 0 and contamination_map > 0:
                    classification_pec = float(classification )/float(contamination_unmap)
                else:
                    classification_pec = 0 
            else: 
                classification = int(merge_reads) - int(unclassified_count)
                if unclassified_count > 0 and host_map > 0:
                    classification_pec = float(classification)/float(merge_reads)
                else:
                    classification_pec = 0 
            if args.contamination or args.neg_control:

                with open ('%s/results/%s.csv' % (results_dir, args.run), 'a') as f:
                    writer = csv.writer(f, delimiter = ',')
                    writer.writerow ([sample, raw_reads, host_map, host_pec,  host_unmap, merge_reads, merge_pec, contamination_map, contamination_unmap, contamination_pec, classification, classification_pec])

            else:
                with open ('%s/results/%s.csv' % (results_dir, args.run), 'a') as f:
                    writer = csv.writer(f, delimiter = ',')
                    writer.writerow ([sample, raw_reads, host_map, host_pec,  host_unmap, merge_reads, merge_pec, classification, classification_pec])
    if args.contamination or args.neg_control:
        for line in fileinput.input(files=['%s/results/%s.csv' % (results_dir, args.run)] , inplace = True):
            if fileinput.isfirstline():
                print 'sample, raw_reads, host_map, host_pec,  host_unmap, merge_reads, merge_pec, contamination_map, contamination_unmap, contamination_pec, classification, classification_pec'
            print line,
    else:
        for line in fileinput.input(files=['%s/results/%s.csv' % (results_dir, args.run)] , inplace = True):
            if fileinput.isfirstline():
                print 'sample, raw_reads, host_map, host_pec,  host_unmap, merge_reads, merge_pec, classification, classification_pec'
            print line,

    print ('information document created, information of reads can be found in %s/results/%s.csv ' % (args.output_dir, args.run))

def predict_genome_cov(args, directory):
# dictionary for genome size
    d_genomesize = defaultdict(list)
    for line in open ('%s/resources/assembly_summary_refseq_less.txt' % (directory)).readlines():
        try:
            line = line.strip().split("\t")
            key, value = line[1], line[0:]
            d_genomesize[key] += value
        except IndexError:
            key, value = 'null', 'null'
            d_genomesize[key] +=value

    d_read_length = defaultdict(str)
    for line in open('%s/multiqc_data/multiqc_general_stats.txt' % (args.output_dir)).readlines():
        line = line.strip().split('\t')
        key,value = line[0], line[3]
        d_read_length[key] += value


# output taxon Id for each sample
    for fastq in os.listdir(args.output_dir):
        if fastq.endswith('.fastq'):
            fastq_split = fastq.split('.')
            sample = fastq_split[0]

            with open ("%s/%s/%s_taxon_ID.csv" %(args.output_dir, sample, sample)  ,  'w' )as f:
                for line in open ("%s/%s/%s_clark_abundance.csv" % (args.output_dir, sample, sample)):
                    line = line.strip()
                    line = line.split("\t")
                    line[0] = line[0].strip()
                    result = line[0], line[2] ,  d_genomesize[line[0]]
                    writer= csv.writer(f, delimiter=',')
                    writer.writerows ([result])

# write clean taxon ID file
    for fastq in os.listdir(args.output_dir):
        if fastq.endswith('.fastq'):
            fastq_split = fastq.split('.')
            sample = fastq_split[0]
            with open("%s/%s/%s_taxon_ID_output.csv" %(args.output_dir,  sample, sample) ,"w") as output_file:
                for line in  open("%s/%s/%s_taxon_ID.csv" %(args.output_dir, sample, sample) ,"rb") :
                    line= line.strip()
                    line = line.split(",")
                    line[2]= line[2].replace("'", "")
                    line[2] = line[2].replace('"', '')
                    line[2]= line[2].replace('[', '')
                    line[2] = line[2].replace(']', '')
                    try:
                        name = line[4]                
                    except IndexError:
                        name = 'error'
                    try:
                        line[6] = line[6].replace("'", "")
                        size = line[6].replace("]", "")
                        size = size.replace('"', '')
        #                print ('size %s' %(size))
                    except IndexError:
                        size= 'null'
                    output = line[0], line[1], size , name

                    writer= csv.writer (output_file , delimiter = ",")
                    writer.writerows ([output])

            os.system ('sort -nr -k2 %s/%s/%s_taxon_ID_output.csv > %s/%s/%s_taxon_ID_outputs.csv' %(args.output_dir,  sample, sample, args.output_dir, sample, sample)) 


    for fastq in os.listdir(args.output_dir):
        if fastq.endswith('.fastq'):
            fastq_split = fastq.split('.')
            sample = fastq_split[0]
            with open("%s/%s/%s_taxon_ID_output_calc.csv" %(args.output_dir, sample, sample) ,"w") as output_file:
                av_length = d_read_length[sample]
                try:
                    av_lenth= float(av_length)
                    print ('av_len = %s ' % (av_length)) 
                    for line in open ("%s/%s/%s_taxon_ID_outputs.csv" %(args.output_dir, sample, sample)).readlines():
                        line= line.strip()
                        line= line.split(",")
                        name=line[3]
                        print name
                        no_reads=line[1]
                        no_read=float(no_reads)
                        print no_reads
                        print no_read
                        try:
                            genome_size_column =line[2]
                            genome_size = float(genome_size_column)
                        except ValueError:
                            genome_size= 'null'
                            print genome_size
                        try:
                            pec_cov= ((no_read*av_lenth)/genome_size)
                        except TypeError:
                            pec_cov = 'null'
                            print pec_cov
                        taxon_ID = line[0]
                        organism = line[3]    
                        writer= csv.writer (output_file , delimiter = ",")
                        writer.writerow ([sample, no_read, genome_size, av_lenth, pec_cov,taxon_ID, organism] )
                except ValueError:                                                                                                                                                                              
                    pass        

    for fastq in os.listdir(args.output_dir):
        if fastq.endswith('.fastq'):
            fastq_split = fastq.split('.')
            sample = fastq_split[0]
            with open ('%s/%s_genome_cov.csv' % (args.output_dir , args.run), 'a') as output_file:
                for line in open ("%s/%s/%s_taxon_ID_output_calc.csv" %(args.output_dir, sample, sample) ).readlines():
                    line= line.strip()
                    line= line.split(",")
                    try: 
                        if float(line[4] )>= 0.001:
                    
                        
                            writer= csv.writer (output_file , delimiter = ",")
                            writer.writerow ([sample, line[6], line[1], line[2], line[3], line[4], line[5]])
                    except ValueError:
                        pass
    os.system('sort %s/%s_genome_cov.csv > %s/%s_genome_cov2.csv' % (args.output_dir, args.run, args.output_dir, args.run))
    os.system('mv %s/%s_genome_cov2.csv %s/%s_genome_cov.csv' % (args.output_dir, args.run, args.output_dir, args.run))

    if args.taxon_contaminants:
        for taxon in open('%s' % (args.taxon_contaminants)):
            taxon = taxon.strip()
            os.system( 'grep %s/%s_genome_cov.csv -v %s > %s/%s_genome_cov2.csv' % (args.output_dir, args.run, taxon, args.output_dir, args.run))
            os.system('mv %s/%s_genome_cov2.csv %s/%s_genome_cov.csv'% (args.output_dir, args.run, args.output_dir, args.run))

    #print ('estimation of genome coverage for each taxon identified in each samples complete, information can be found in %s/%s_genome_cov.csv' % (args.output_dir, args.run))

def auto_assemble(args, directory):
    count = 0
    for line in open ('%s/%s_genome_cov.csv' % (args.output_dir, args.run)).xreadlines(): count +=1
    if count > 1:
        out_dir = os.path.abspath('%s' % (args.output_dir))
        d_read_in_ftp = defaultdict(str)

        for line in open("%s/resources/assembly_summary_refseq.txt" % (directory)).readlines():
            try:
                line = line.strip()
                line = line.split("\t")
                
                key, value = line[5], line[19]
                d_read_in_ftp[key] += value
            except IndexError:
                key, value = 'null', 'null'
                d_read_in_ftp[key] += value    

        with open ('%s/%s_sample_ref.csv' % (args.output_dir, args.run), 'w') as f:

            for line in open ("%s/%s_genome_cov.csv" % (args.output_dir, args.run), 'rb'):
                line = line.strip().split(',')
                line[5] = float(line[5])
               
                line[3] = line[3].strip()
                line[3] = str(line[3])
                genome_size = line[3]
           
                line[6]=line[6].strip()
                line[6]= str(line[6])
                result = d_read_in_ftp[line[6]]
                result_split = result.split('/')
            
                try: 
                    result_ftp = ['ftp://' + result_split[-8]+'/' + result_split[-7]    +'/' + result_split[-6]+ '/' + result_split[-5]+'/' + result_split[-4]+'/'+ result_split[-3]+'/' + result_split[-2]+'/' + result_split[-1]+'/' + result_split[-1] + '_genomic.fna.gz']
                except IndexError:
                    print 'error'
                result_ftp = str(result_ftp)       
                result_ftp = result_ftp.replace('[', '')
                result_ftp = result_ftp.replace(']', '')    
                result_ftp  = result_ftp.replace("'", "")
                result_ftp_unzip = [result_split[-1] + '_genomic.fna']
                result_ftp_unzip = str(result_ftp_unzip)
                result_ftp_unzip = result_ftp_unzip.replace('[', '')
                result_ftp_unzip = result_ftp_unzip.replace(']', '')
                result_ftp_unzip = result_ftp_unzip.replace("'", "")
                result_filename = [result_split[-1] + '_genomic.fna.gz']
                result_filename = str(result_filename)
                result_filename = result_filename.replace('[' , '')
                result_filename = result_filename.replace(']' , '')
                result_filename = result_filename.replace("'", "")
                
                name = line[1].replace("'", "")
                name = name.replace(" ", "_")
                name = name.replace("[", "")
                name = name.replace("]", "")
                
                try:
                    os.mkdir('%s/refs' % (directory))
                except OSError:
                    pass
                    
                #print ('%s is ref' % (result_filename))                      
                #if os.path.exists('%s/refs/%s' % (directory, result_filename)):
                 #   print 'exists'
                  #  os.system ('gunzip -f %s/refs/%s' % (directory, result_filename ))
                writer= csv.writer (f, delimiter = ',')
                writer.writerow ([line[0],name, genome_size, line[2], line[5], line[6] , result_ftp_unzip , result_ftp])
                    
                #if not os.path.exists('%s/refs/%s' % (directory, result_filename)):
                 #   print line[0]  
                  #  print ('downloading %s' % (result_filename))  
                    #os.system ('wget -P %s/refs %s ' % (directory, result_ftp))
                   # os.system ('gunzip -f %s/refs/%s' % (directory, result_filename ))
                    
            
                    #writer= csv.writer (f, delimiter = ',')
                    #writer.writerow ([line[0],name, genome_size, line[2], line[5], line[6] , result_ftp_unzip , result_ftp])

        with open ('%s/%s_bwa_index.csv' % (args.output_dir, args.run) , 'a') as bf:
            for line in open ('%s/%s_sample_ref.csv' % (args.output_dir, args.run )).readlines():
                line= line.strip()
                line = line.split(",")
                ref = line[6]
                ref_path = ('%s/refs/%s' %  (directory, ref) 
                #ref_path =  ref_path.replace('[' , '')
                #ref_path =  ref_path.replace(']' , '')
                #ref_path =  ref_path.replace("'", "")

                ref_ftp = line[7]
                if os.path.exists ('%s/refs/%s.pac' % (directory, ref)):
                    print 'reference indexed'
                else: 
                    writer = csv.writer (bf, delimiter = ',')
                    writer.writerow ([ref_path, ref_ftp])
                    print ref_path
        for line in open ('%s/%s_bwa_index.csv' % (args.output_dir, args.run) ).readlines():
            line = line.strip().split(',')
            ref_filename = line[1]
            #ref_filename = ref_filename.replace('[' , '')
            #ref_filename = ref_filename.replace(']' , '')
            #ref_filename = ref_filename.replace("'", "")
            ref_gunzip = line[0]
            #ref_gunzip = ref_gunzip.replace('[' , '')
            #ref_gunzip = ref_gunzip.replace(']' , '')
            #ref_gunzip = ref_gunzip.replace("'", "")
            ref_gunzip = ref_gunzip + '.gz'

           
            os.system ('wget -P %s/refs %s ' % (directory, ref_filename))
            os.system ('gunzip %s' % (ref_gunzip)) 
# this is looping twice .......
        

        os.system ('nextflow run %s/nextflow_scripts/bwa_index.nf --refdir %s/refs' % (directory , directory))
        try:
            os.mkdir('%s/reference_mapping' % (args.output_dir))
        except OSError:
            pass
        os.system ('nextflow run %s/nextflow_scripts/bwa_mem.nf --threads %s --list  %s/%s_sample_ref.csv --path %s/paired/ --refdir %s/refs/ --dataDir %s/reference_mapping/' % (directory, args.threads, out_dir, args.run, out_dir, directory, out_dir))

        #os.system ('gzip %s/refs/*.fna' % (directory)) 
        
        with open ('%s/%s_reference_mapping.csv' % (args.output_dir, args.run) , 'a') as f:
            for each_line in open ('%s/%s_sample_ref.csv' % (args.output_dir, args.run)):
                each_line = each_line.strip().split(',')
                sample = each_line[0]
                ref = each_line[6]
                genome_size = each_line[2].strip()
                genome_size = float (genome_size)
                organism = each_line[1]
                clark_reads = each_line[3]
                clark_pec = each_line[4]
                taxon_ID = each_line[5]


                read_count = 0
                for line in open ('%s/paired/%s_R1.fastq' % (args.output_dir, sample)).xreadlines( ): read_count +=1
                reads = read_count/4

                for line in open ('%s/reference_mapping/%s%s.txt' % (args.output_dir, sample, organism)):
                    line = line.strip()
                    genome_cov = float(line)

                for line in open ('%s/reference_mapping/%s%s_mapped.txt' % (args.output_dir, sample, organism)):
                    line = line.strip()
                    mapped_reads = float(line)   

                pec_genome_cov = genome_cov/genome_size

                pec_mapped_reads = mapped_reads/reads

                writer = csv.writer(f, delimiter = ',')
                writer.writerow([sample, reads, organism, genome_size, genome_cov, pec_genome_cov, mapped_reads, pec_mapped_reads, taxon_ID, ref,clark_reads, clark_pec])
        os.system ('sort %s/%s_reference_mapping.csv > %s/%s_reference_mappings.csv' % (args.output_dir, args.run, args.output_dir, args.run) )
        os.system ('mv %s/%s_reference_mappings.csv  %s/%s_reference_mapping.csv' % (args.output_dir, args.run, args.output_dir, args.run) )

        for line in fileinput.input(files=['%s/%s_reference_mapping.csv' % (args.output_dir, args.run)] , inplace = True):
            if fileinput.isfirstline():
                print 'sample, reads, organism,genome_size,genome_cov,pec_genome_cov,mapped_reads,pec_mapped_reads,taxon_ID,ref,clark_reads, clark_pec'
            print line,

        print ('automatic assembly of samples against identified refernce completed, information can be found in %s/%s_reference_mapping.csv' % (args.output_dir, args.run))
    else:
        print ('no significant genoomes classified, no refernce assembly performed')
    os.system ('gzip -f %s/refs/*.fna' %(directory))    
def taxon_abundance(args):
    d_taxon_reads = {}

    with open('%s/%s_genome_cov.csv' % (args.output_dir, args.run), 'r') as data_file:
        for row in data_file:
            row= row.strip().split(',')
            taxonID = row[6]
            sample = row[0]
            reads = float(row[2])

            if taxonID in d_taxon_reads:
                if sample in d_taxon_reads[taxonID]:
                    print 'error in d_taxon_reads'
                else:
                    d_taxon_reads[taxonID]
                    d_taxon_reads[taxonID][sample] = reads

            else:
                d_taxon_reads[taxonID] = {}
                d_taxon_reads[taxonID][sample] = reads

    print d_taxon_reads
    d_taxon_genome_cov = {}
    count = 0
    if os.path.isfile('%s/%s_reference_mapping.csv' % (args.output_dir, args.run)):
        for line in open ('%s/%s_reference_mapping.csv' % (args.output_dir, args.run)).xreadlines(): count +=1
        if count > 1:

            with open ('%s/%s_reference_mapping.csv' % (args.output_dir, args.run), 'r') as data:
                next(data)
                for line in data:
                    line = line.strip().split(',')
                    taxon = line[8]
                    sample = line[0]
                    try:
                        genome_cov = float(line[4])

                    except ValueError:
                        pass
                    if taxon in d_taxon_genome_cov:
                        if sample in d_taxon_genome_cov[taxon]:
                            print 'error in d_taxon_genome_cov'
                        else:
                            d_taxon_genome_cov[taxon]
                            d_taxon_genome_cov[taxon][sample] = genome_cov

                    else:
                        d_taxon_genome_cov[taxon] = {}
                        d_taxon_genome_cov[taxon][sample] = genome_cov
            print d_taxon_genome_cov
        else:
            pass
    else:
        pass
    d_taxon_organism = {}            
    with open ('%s/%s_genome_cov.csv' % (args.output_dir, args.run), 'r') as data:
        for line in data:
            line = line.strip().split(',')
            taxon = line[6]
            
            organism = line[1]
            
            key, value = taxon, organism
            d_taxon_organism[key] =value
    

    with open ('%s/%s_taxonID_info.csv' % (args.output_dir, args.run), 'w') as f:
        for key in d_taxon_reads.iterkeys():
            i=0
            for item in d_taxon_reads[key]:
                if len(item) > 0:
                    i=i +1

            L = sorted(d_taxon_reads[key].values())
            m = i - 1
            median= (L[i/2] + L[m/2]) / 2.0        
            taxonID_number =i
            Taxon_ID= key
            sample_reads =d_taxon_reads[key]
            number_of_reads = d_taxon_reads[key].values()      
            total_reads= sum(number_of_reads)
            average_reads =total_reads/i
            highest_reads= max(number_of_reads) 
            for sample in d_taxon_reads[key]:
                if d_taxon_reads[key][sample] == highest_reads:
                    print key, highest_reads, sample, d_taxon_reads[key][sample]
                    max_sample = sample
            print max_sample
            lowest_reads= min(number_of_reads)
        
            try:
                genome_cov = d_taxon_genome_cov[key].values()
                max_cov = max(genome_cov)
                min_cov = min(genome_cov)
                total_cov = sum(genome_cov)
                n = 0 
                for thing in d_taxon_genome_cov[key]:
                    if len(thing) > 0:
                        n = n+1

                av_cov = total_cov/n
                X = sorted (genome_cov)
                p = n -1
                if max_cov > 0:
                    median_cov = (X[n/2] + X[p/2]) /2 
                else:
                    median_cov = 0
                for covs in d_taxon_genome_cov[key]:
                    if d_taxon_genome_cov[key][covs] == max_cov:
                        max_cov_sample = covs
            except KeyError:
                genome_cov = 'NA'
                av_cov = 'NA' 
                min_cov = 'NA'
                total_cov = 'NA'
                max_cov = 'NA'
                median_cov = 'NA'
                max_cov_sample = 'NA'
            organism_name = d_taxon_organism[key]
            writer= csv.writer(f, delimiter=',')
            writer.writerow ([Taxon_ID,  organism_name, taxonID_number, average_reads, median, lowest_reads, highest_reads, max_sample,  av_cov, median_cov, min_cov, max_cov, max_cov_sample])

        for line in fileinput.input(files=['%s/%s_taxonID_info.csv' % (args.output_dir, args.run)] , inplace = True):
            if fileinput.isfirstline():
                print 'Taxon_ID,  organism_name, taxonID_number, average_reads, median, lowest_reads, highest_reads, max_sample,  av_cov, median_cov, min_cov, max_cov, max_cov_sample'
            print line,
    os.system('cp %s/*.csv %s/results/' %(args.output_dir, args.output_dir))
    os.system('cp %s/*/*.html %s/results' % (args.output_dir, args.output_dir))   

    print ('report on frequency at which taxon IDs were identified in the run complete, see %s/%s_taxonID_info.csv ' % (args.output_dir, args.run))    

def remove_work_dir ():
    os.system ('rm -r work')

main()