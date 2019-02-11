import glob
import concurrent.futures
import os
import fileinput
import csv 


def header(each_file):

	with  open (each_file, 'rU') as f:
        name = each_file 
        for line in f.readlines():
                    #print line  
            if line.startswith('>'):
                    #print line           
                    
                split_file = line.split('|')
                    #print split_file[2]
                old = line
                print old
                acc = split_file[2]
                ID = split_file[3]
                header = '>' +  acc 
            with open (each_file, 'U') as fi:    
                newText = fi.read()
                while ('%s' %( old)) in newText:
                    newText=newText.replace ('%s' %(old) , '%s \n' % (header))

            with  open (each_file, 'w') as fi:
                fi.write (newText)
            os.system ("sed -i '/^$/d' %s" % (each_file) )

with concurrent.futures.ProcessPoolExecutor() as executor:
    each_file = glob.glob('*.fasta')

    executor.map(header, each_file)