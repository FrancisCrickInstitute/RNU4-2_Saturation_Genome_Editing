#run_seqprep_210729.py

#compatible with sample names that still contain the "_S" from demultiplexing with bcl2fastq on our own
#to be run in the folder for the experiment with all the reads split by sample as fastq.gz files
#generates a shell script to run
#requires no underscores to be in sample names at this point
#outputs the stats from each seqprep call to a file called seqprep_stats.txt

'''Adapter for PU1L , PU1R:
-A GGTTTGGAGCGAGATTGATAAAGT #if PU1R is used (as seen in R1):  
-B CTGAGCTCTCTCACAGCCATTTAG #if PU1L is used (as seen in R2): 

For ASF Nextseq adapters:
-A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC (as seen in R1)
-B CTGTCTCTTATACACATCTGACGCTGCCGACGA (as seen in R2)


'''

#usage python run_seqprep_VHL_pipeline.py [adapter1] [adapter2]
import os
import sys
import subprocess

adapter1 = sys.argv[1]
adapter2 = sys.argv[2]
working_dir = os.getcwd()
command_file_name = working_dir+"/run_seqprep.sh"
command_file = open(command_file_name, 'w')
command_file.write('module load SeqPrep\n')
for i in os.listdir(os.getcwd()):
    if i.startswith("Undetermined"):
        pass
    elif i.endswith("R1_001.fastq.gz"): 
        index_of_first_under = i.find('_')
        index_of_second_under = i.find('_',index_of_first_under+1)
        index_of_third_under = i.find('_',index_of_second_under+1)
        index_of_sample_number = i.find('_S')+1
        sample_name = i[:index_of_second_under]
        sample_name_dash = i[:index_of_first_under]+'-'+i[index_of_first_under+1:index_of_second_under]
        file_name = i[:index_of_third_under]
        # here are adapters for F and R seq primers
        command_file.write(r'SeqPrep -f '+file_name+r'_R1_001.fastq.gz -r '+file_name+r'_R2_001.fastq.gz -1 Seqprep/R1/'+sample_name_dash+r'.R1.fastq.gz -2 Seqprep/R2/'+sample_name_dash+r'.R2.fastq.gz -A '+adapter1+' -B '+adapter2+' -M 0.1 -s Seqprep/merged/'+sample_name_dash+'.merged.fastq.gz -m 0.001 -q 20 -o 15 &\n')
    else:
        pass
command_file.write("wait\n")
command_file.close()



