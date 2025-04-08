#run_remove_n_bases.py
#modified to now require a no_Ns directory in place to put merged files (without Ns)
#output will be unzipped fastq's
import os
import sys
import subprocess
#hard-coded directories
#will output a script to be run in shell

working_dir = os.getcwd()
shell_file = open('run_remove_n_bases.sh', 'w')
shell_file.write('cd '+working_dir+'\n')

for i in os.listdir(working_dir):
    if i.endswith(".fastq.gz"): #will only work on gzipped fastq's
        index_of_first_dot = i.find('.')
        index_of_first_under = i.find('_')
        index_of_extension = i.find('.fastq.gz')
        sample_name = i[:index_of_first_dot]
        sample_lib = i[:index_of_first_under]
        before_extension_name = i[:index_of_extension]
        shell_file.write("python ~/home/users/findlag/bin/remove_n_bases.py "+i+' '+'no_Ns/'+before_extension_name+'.fastq &\n')
    else:
        pass

shell_file.write('wait\n')
shell_file.close()

#script format:  python