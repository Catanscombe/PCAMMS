import argparse
import os
import csv
import sys
import fileinput

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument( 'input_dir',  help = 'directory containing input files')
    parser.add_argument('output_dir', help ='output directory')
    parser.add_argument( 'run',  help = 'run description')
    parser.add_argument ('-nl' , '--neg_list' , help = 'path to csv containing information on negative control samples, see example for format')
    parser.add_argument('-r' , '--rename' , choices = ['yes' , 'no'], default = 'yes', help = 'option to rename files when in results directory, based on sample number and run name, default is yes')
    parser.add_argument('-id' , '--sample_ID', help = 'path to csv containing sample number and sample ID, will rename samples based on sample number and sample ID, see example for format. If the run contains one or more negative control, please submit seperately using --neg_list command' )
    parser.add_argument('-s' , '--sample_list' , help = 'provide a list of samples to process, default is to process all files in input directory')
    parser.add_argument('-n1', '--neg_sample_R1',  help = 'path to negative control file R1')
    parser.add_argument('-n2', '--neg_sample_R2', help = 'path to negative control file R2')
    parser.add_argument ('-t' , '--threads' , default = '4', help = 'number of threads available, default is 4')
    parser.add_argument( '-rh', '--remove_host', choices = ['yes', 'no'], default = 'yes' , help = 'option remove host reads at the beginning  of analysis, default yes. Default host is Human, if another host is needed, place the reference   fasta in host directory. See README for more detail' )

    args = parser.parse_args()
    
    directory = os.path.dirname(os.path.abspath(sys.argv[0]))

    id_dir (args, directory)
    rename(args)
    remove_host(args, directory)
    qc_reads(args, directory)
    neg_control_class (args, directory)
    remove_work_dir ()

def id_dir( args, directory):
## check if input and output directorys exists, if optional files are given, checks the file path,  makes sample list
    #check input dir exists
    counter =0

    if not  os.path.exists(args.output_dir):
        try:
            os.mkdir('%s' % (args.output_dir))
        except OSError :
            print ('Error! cannot create results directory %s, please check path' % (args.output_dir))  
            counter = counter + 1
        else:
            print ('%s created, results will be written to this directory' % (args.output_dir))
    else:
        print ('Output directory %s already exisits results will be written to this directory' % (args.output_dir))

    #check output dir
    if  os.path.isdir(args.input_dir):
        print ('input_dir ok')
        if args.sample_list:
            os.system('cp %s %s/sample_list_%s.csv' % (args.sample_list, args.output_dir, args.run))
        else:
            input_list =(f for f in os.listdir(args.input_dir) if f.endswith(('gz' , 'fq' ,'fastq')))
            with open ('%s/sample_list_%s.csv'% (args.output_dir, args.run), 'wb') as out:
                writer = csv.writer(out)
                list_input = zip(input_list)
                writer.writerows(list_input)
    else:    
        print ('Error!  %s is not a valid directory, please provide valid directory for input files' % (args.input_dir))
        counter = counter + 1
    
    #check file containing negative control list exists (if given)   
    if args.neg_list:
        if os.path.isfile(args.neg_list):
            print 'negative control list given'
        else: 
            print 'Error! file containing negative control cannot be found, please check path'
            counter = counter + 1

    #check negative control files (if given)
    if args.neg_sample_R1:
        if os.path.isfile(args.neg_sample_R1):
            print 'negative control R1 ok'
            if os.path.isfile(args.neg_sample_R2):
                print 'negative control ok'
            else:
                print 'problem with negative control R2'
                counter = counter +1 

        else:
            print 'Error! negative control files cannot be found'
            counter = counter +1 

    #will files be renamed?    
    if args.rename == 'yes':
        print 'sample names will be changed in results directory'

    #check sample list (if given)
    if args.sample_list:
        if os.path.isfile(args.sample_list):
            print 'sample list given, will process samples on list within input directory'
        else:
            print 'Error! cannot find sample list, please check path'
            counter = counter + 1
    # check sample ID list if given)
    if args.sample_ID:
        if os.path.isfile(args.sample_ID):
            print 'sample ID list provided, samples will be renamed according to the file'
        else:
            print 'Error! sample list cannot be found, please check path'
            counter = counter + 1
    if counter >= 1:
        print ('problems with %s  input file or directory' %(counter))
        sys.exit()
    print ('this script is runing from %s' % (directory))
   

def rename(args):
## copy paired end reads to paired directory in output directory and rename where appropriate

    # make a dir in results to house paried end reads    
    try:
        os.mkdir('%s/paired' % (args.output_dir))
    except OSError:
        pass
    
    paired_dir = (args.output_dir + "/paired")

    # move and rename negative control if single sample given (-n1/2 or --neg_sample_1/2 command), rewrite sample list to remove this sample from further renaming. 
    
    if args.neg_sample_R1:
        neg_directory = os.path.dirname(os.path.abspath(args.neg_sample_R1))
        print ('negative controls will be processed in %s' %  (neg_directory))
        
        neg_sample_R1 = args.neg_sample_R1
        neg_sample_R1_split = neg_sample_R1.split('/')
        neg_sample_R1_ID = neg_sample_R1_split[-1]
        
        neg_sample_R2 = args.neg_sample_R2
        neg_sample_R2_split = neg_sample_R2.split('/')
        neg_sample_R2_ID = neg_sample_R2_split[-1]
        
        neg_controls = ([neg_sample_R1_ID, neg_sample_R2_ID])
        with open ('%s/sample_list_%s.csv' % (args.output_dir, args.run)) as oldfile, open ('%s/sample_list_%s2.csv' % (args.output_dir, args.run), 'w') as newfile:
            for line in oldfile:
                if not any( neg_control in line for neg_control in neg_controls ):
                    newfile.write(line)
        
        if neg_sample_R1_ID.endswith ('gz'):
            os.system ('cp %s %s/neg_%s_R1.fastq.gz' % (neg_sample_R1, paired_dir , args.run))
            os.system ('gunzip  %s/neg_%s_R1.fastq.gz' % (paired_dir , args.run))
        else: 
            os.system ('cp %s %s/neg_%s_R1.fastq' % (neg_sample_R1, paired_dir , args.run))

        if neg_sample_R2_ID.endswith ('gz'):
            os.system ('cp %s %s/neg_%s_R2.fastq.gz' % (neg_sample_R1, paired_dir , args.run))
            os.system ('gunzip  %s/neg_%s_R2.fastq.gz' % (paired_dir , args.run))
        else:
            os.system ('cp %s %s/neg_%s_R2.fastq' % (neg_sample_R2, paired_dir , args.run))
        
        os.system ('mv %s/sample_list_%s2.csv %s/sample_list_%s.csv' % (args.output_dir, args.run, args.output_dir, args.run))   
        print ('negative controls copied  to %s and renamed' % (paired_dir))     
        
    # rename based on input IDs from user (-id --sample_id)
        
    if args.sample_ID: 
        with open ('%s' % (args.sample_ID), 'rb') as f:
            next(f)
            for line in f:
                line= line.strip()
                line_split =line.split(',')
                print line_split
                sample_ID = line_split[0]
                sample_number= line_split[1]
                for each_file in os.listdir('%s' % (args.input_dir)):
                    if each_file.startswith('%s_' % (sample_number)):
                        if '_R1' in each_file:
                            if each_file.endswith('.gz'):
                                os.system ('cp %s/%s %s/%s_%s_R1.fastq.gz' % (args.input_dir, each_file, paired_dir, sample_ID, args.run))
                                os.system ('gunzip %s/%s_%s_R1.fastq.gz' % (paired_dir, sample_ID, args.run) )
                            else:
                                os.system ('cp %s/%s %s/%s_%s_R1.fastq' % (args.input_dir, each_file, paired_dir, sample_ID, args.run))

                        elif ('%s_' % (sample_number)) and '_R2' in each_file:
                            if each_file.endswith('.gz'):
                                os.system ('cp %s/%s %s/%s_%s_R2.fastq.gz' % (args.input_dir, each_file, paired_dir, sample_ID, args.run))
                                os.system ('gunzip %s/%s_%s_R2.fastq.gz' % (paired_dir, sample_ID, args.run) )
                            else:
                                os.system ('cp %s/%s %s/%s_%s_R2.fastq' % (args.input_dir, each_file, paired_dir, sample_ID, args.run)) 
        
        print ('samples coped to %s and renamed based on sample ID input from %s' % (paired_dir, args.sample_ID))

    # or option to automatically rename files to number_runid
    elif args.rename == 'yes':
        for line in open ('%s/sample_list_%s.csv' % (args.output_dir, args.run)):
            line = line.strip()
            line_split= line.split('_')
            if '_R1' in line:
                if line.endswith('.gz'):
                    os.system ('cp %s/%s %s/%s_%s_R1.fastq.gz' % (args.input_dir,line, paired_dir, line_split[0], args.run))
                    os.system ('gunzip %s/%s_%s_R1.fastq.gz' % (paired_dir, line_split[0], args.run) )
                else: 
                    os.system ('cp %s/%s %s/%s_%s_R1.fastq' % (args.input_dir,line, paired_dir, line_split[0], args.run))
            if '_R2' in line:
                if line.endswith('.gz'):
                    os.system ('cp %s/%s %s/%s_%s_R2.fastq.gz' % (args.input_dir,line, paired_dir, line_split[0], args.run))
                    os.system ('gunzip %s/%s_%s_R2.fastq.gz' % (paired_dir, line_split[0], args.run) )
                else: 
                    os.system ('cp %s/%s %s/%s_%s_R2.fastq' % (args.input_dir,line, paired_dir, line_split[0], args.run))
        
        print ('samples copied to %s and renamed based on file number and run ID' % (paired_dir))
    
    else:
        for line in open ('%s/sample_list_%s.csv' % (args.output_dir, args.run)):
            line = line.strip()
            line_split= line.split('.')
            if line.endswith('.gz'):
                os.system('cp %s/%s %s/%s.fastq.gz' % (args.input_dir,line, paired_dir, line_split[0]))
                os.system('cp %s/%s %s/%s.fastq.gz' % (args.input_dir,line, paired_dir, line_split[0]))
            else:
                os.system('cp %s/%s %s/%s.fastq' % (args.input_dir,line, paired_dir, line_split[0]))

        print ('samples copied to %s'  % (paired_dir)) 

    # rename negative controls when multiple are given
    if args.neg_list:
        if not args.sample_ID:
            with open ( '%s' % (args.neg_list), 'rb') as f:
                next(f)
                for line in f:
                    split_line = line.strip().split(',')
                    Id = split_line[0]
                    number = split_line[1]
                    print Id
                    print number
                    os.system ('mv %s/paired/%s_%s_R1.fastq %s/paired/%s_%s_R1.fastq' % (args.output_dir , number, args.run, args.output_dir, Id, args.run))
                    os.system ('mv %s/paired/%s_%s_R2.fastq %s/paired/%s_%s_R2.fastq' % (args.output_dir , number, args.run, args.output_dir, Id, args.run))
                print ('negative controls copied to  %s and renamed base on %s' % (paired_dir, args.neg_list))
    
    if args.neg_list and args.sample_ID:
        with open ( '%s' % (args.neg_list), 'rb') as f:
            next(f)
            for line in f:
                split_line = line.strip().split(',')
                Id = split_line[0]
                number = split_line[1]

                print Id
                print number
                for each_file in os.listdir('%s' % (args.input_dir)):
                    if each_file.startswith('%s_' % (number)):
                        if each_file.endswith('.gz'):
                            os.system('gunzip %s' % (each_file))
                        if '_R1' in each_file:
                            os.system ('cp %s/%s %s/paired/%s_%s_R1.fastq' % (args.input_dir , each_file, args.output_dir, Id, args.run))
                        elif '_R2' in each_file:
                            os.system ('cp %s/%s %s/paired/%s_%s_R2.fastq' % (args.input_dir , each_file, args.output_dir, Id, args.run))
                print ('negative controls copied to %s and renamed base on ' % (paired_dir, args.neg_list))
        
    
def remove_host(args, directory):   
## remove reads which map to host genome, default is on and human, host can be changed, see readme file
    
    if args.remove_host == 'yes':
        # make directory for qc result storage
        try:
            os.mkdir('%s/qc_results' % (args.output_dir))
        except OSError:
            pass
            print directory 

        # refernce map to host genome (bwa mem) then use samtools to get mapping details eg number of reads mapped, create a new fastq file from unmapped reads
        for fasta in os.listdir('%s/host'% (directory)):
            if fasta.endswith('.fasta.gz'):
                print fasta
                host_path = os.path.dirname(os.path.abspath('%s/host/%s'% (directory , fasta)))
                print host_path
                try: 
                    os.mkdir('%s/paired/map_host' % (args.output_dir))
                except OSError:
                    pass
                if args.threads < 4:
                    os.system ('nextflow run %s/nextflow_scripts/host_map.nf --host_path "%s/%s" --outdir "%s/paired/" --threads "%s"' % (directory, host_path, fasta, args.output_dir, args.threads))
                else: 
                    os.system ('nextflow run %s/nextflow_scripts/host_map.nf --host_path "%s/%s"  --outdir "%s/paired/" ' % (directory, host_path, fasta, args.output_dir))

        # make a csv with details of the host map process
        with open ('%s/qc_results/%s_host_map_info.csv' % (args.output_dir, args.run), 'a') as f:
            writer = csv.writer (f, delimiter = ',')
            for each_file in os.listdir ('%s/paired' % (args.output_dir)):
                if each_file.endswith('_R1.fastq'):
                    each_file_split = each_file.split('.')
                    sample_R1 = each_file_split[0]
                    sample = sample_R1.replace('_R1', '')
                    raw_count = 0
                    for line in open ('%s/paired/%s' % (args.output_dir, each_file)).xreadlines( ): raw_count +=1
                    host_unmap_count = 0
                    for line in open ('%s/paired/map_host/%s' % (args.output_dir,each_file)).xreadlines(): host_unmap_count +=1
                    raw_reads = int(raw_count/4)
                    host_unmap_reads = host_unmap_count/4
                    host_map_reads = raw_reads - host_unmap_reads
                    if host_map_reads > 0 and raw_reads > 0:
                        pec_host_map = float(host_map_reads)/float(raw_reads)
                    else:
                        pec_host_map = 0
                    writer.writerow([sample, raw_reads, host_map_reads , pec_host_map, host_unmap_reads])

        # sort the file by sample and then add headings to columns
        os.system('sort %s/qc_results/%s_host_map_info.csv > %s/qc_results/%s_host_map_info2.csv' % (args.output_dir, args.run, args.output_dir, args.run))
        os.system ('mv %s/qc_results/%s_host_map_info2.csv  %s/qc_results/%s_host_map_info.csv' % (args.output_dir, args.run, args.output_dir, args.run))
        for line in fileinput.input(files=['%s/qc_results/%s_host_map_info.csv' % (args.output_dir, args.run)] , inplace = True):
                if fileinput.isfirstline():
                    print 'sample,raw_reads,host_map_reads,pec_host_map,output_reads'
                print line,
        os.system('cp %s/paired/map_host/*.fastq %s/paired/' % (args.output_dir, args.output_dir))
        os.system('rm -r  %s/paired/map_host' % (args.output_dir))

def qc_reads(args, directory):
## merge paried end reads, and perform basic qc on merged reads
    
    
    # perform fastqc on all merged reads, then perform mulitqc
    os.system ('nextflow run %s/nextflow_scripts/qc.nf --outdir "%s/paired/" ' % (directory,  args.output_dir))
    print 'fastqc performed on  reads'
    os.system ('mv %s/paired/fastqc_output %s/fastqc_output' % (args.output_dir, args.output_dir    )) 
    os.system ('multiqc %s/fastqc_output -o %s' %(args.output_dir, args.output_dir))
    print 'multiqc performed on all  reads'

    # create a csv with read information in, including number of raw forward reads (after host removal and how many reads were sucessfully merged. 
    with open ('%s/%s_read_info.csv' % (args.output_dir, args.run), 'a') as f:
        writer = csv.writer (f, delimiter = ',')
        for each_file in os.listdir ('%s/paired' % (args.output_dir)):
            if each_file.endswith('.fastq'):
                each_file_split = each_file.split('_')
                sample = each_file_split[:-1]
                sample ='_'.join(sample)
                print sample
                R1_count = 0
                for line in open ('%s/paired/%s_R1.fastq' % (args.output_dir, sample)).xreadlines(): R1_count +=1
                R2_count = 0
                for line in open ('%s/paired/%s_R2.fastq' % (args.output_dir, sample)).xreadlines(): R2_count +=1
                R1_reads = R1_count/4
                R2_reads = R2_count/4
                
                writer.writerow([sample,  R1_reads, R2_reads])
    os.system ('sort %s/%s_read_info.csv > %s/%s_read_info2.csv' % (args.output_dir, args.run, args.output_dir, args.run))
    os.system (' mv %s/%s_read_info2.csv %s/%s_read_info.csv' % (args.output_dir, args.run, args.output_dir, args.run))          
    #os.system ('gzip %s/paired/*.fastq' % (args.output_dir))
    for line in fileinput.input(files =['%s/%s_read_info.csv' % (args.output_dir, args.run)] , inplace = True):
        if fileinput.isfirstline():
            print 'sample,  R1_reads, R2_reads'
        print line, 

def neg_control_class (args, directory):
## If negative control(s) were given, classify these samples using CLARK-L and produce a Krona plot. 
    
    #get the directory that contains the clark scripts and database
    Clark_dir = ('%s/CLARKSCV1.2.6' % (directory))
    
    try:
        os.mkdir('%s/qc_results' % (args.output_dir))
    except OSError:
        pass
    os.system ('mv %s/*.html %s/qc_results' % (args.output_dir, args.output_dir))
    os.system ('mv %s/%s_read_info.csv %s/qc_results ' % (args.output_dir, args.run, args.output_dir))
    
    #if a single negative control was given, make a negative control directoy and move merged reads into this directory 
    if args.neg_sample_R1: 
        try:
            os.mkdir ('%s/neg_control' % (args.output_dir))
        except OSError:
            pass

        os.system ('mv %s/paired/neg_%s*.fastq %s/neg_control/' % (args.output_dir, args.run, args.output_dir))
        
    #if a list of  negative controls was given, make a negative control directoy and move merged reads into this directory  
    if args.neg_list:
        try:
            os.mkdir('%s/neg_control' % (args.output_dir))
        except OSError:
            pass
        for line in open (args.neg_list):
            split_line = line.strip().split(',')
            Id = split_line[0]
            number = split_line[1]
            os.system ('mv %s/paired/%s_%s*.fastq %s/neg_control/' % (args.output_dir, Id, args.run, args.output_dir))
    # If any negative controls were given classify these reads using clark and make a krona plot
    if args.neg_sample_R1 or args.neg_list:
        for files in os.listdir('%s/neg_control' % (args.output_dir)):
            if files.endswith('_R1.fastq'):
                neg_control_directory = os.path.dirname(os.path.abspath('%s/neg_control/%s'% (args.output_dir,  files)))
                split_files = files.split('_')
                sample = split_files[:-1]
                sample ='_'.join(sample)
                print sample
                cwd = os.getcwd()
                os.chdir ('%s' % (Clark_dir))
                os.system ( './classify_metagenome.sh  --light -P  %s/%s_R1.fastq %s/%s_R2.fastq  -n %s -R %s/%s_clark' % ( neg_control_directory, sample, neg_control_directory, sample, args.threads , neg_control_directory, sample))
                os.system ('./estimate_abundance.sh --krona -F %s/%s_clark.csv -D DIR_DB > %s/%s_clark_abundance.csv' % ( neg_control_directory, sample, neg_control_directory,  sample))
                os.system ('mv results.krn %s/%s_clark_abundance.csv' % ( neg_control_directory,  sample ))
                os.chdir ('%s' % (cwd))
                os.system  ('ktImportTaxonomy -o  %s/%s.html -m 3 %s/%s_clark_abundance.csv' % ( neg_control_directory, sample, neg_control_directory, sample))
                
                #move results to qc_results directory 
                os.system (' cp %s/neg_control/*.html %s/qc_results/' % (args.output_dir, args.output_dir))
        print 'classification performed on negative control(s)'


def remove_work_dir ():
## remove the work directory created by next flow
    os.system ('rm -r work')
main()

