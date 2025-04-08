#run_needle_to_sam_BRCA1_pipeline.py
#will output a script to be run in shell
#to be run in directory with the fastq files (merged, unzipped)
#this re-works naming scheme to switch from dashes to UNDERSCORES
#uses the dot in the .merged.fastq to remove suffix / define sample 'i.e. X17r1-pre.merged.fastq'

#amplicon+replicate+'-'+timepoint.merged.fastq

import os
import sys
import subprocess

working_dir = os.getcwd()
shell_file = open('run_needle_to_sam.sh', 'w')
shell_file.write('module load EMBOSS/6.6.0-foss-2016b\n') #ideally would implement version control here...

ref_folder = sys.argv[1]
#reference comes first, and then the cigar reflects changes from reference to sample
for i in os.listdir(working_dir):
    if i.endswith(".fastq"): 
        index_of_first_dot = i.find('.')
        sample_name = i[:index_of_first_dot]
        index_of_first_dash = i.find('-')
        if 'r' in sample_name: #for i.e. X17r1-pre
            sample_amplicon = sample_name[:i.find('r')]
        else:
            #for i.e. X17-lib or X5-neg
            sample_amplicon = i[:index_of_first_dash]
        sample_ref = ref_folder+'/'+sample_amplicon+'.fa'
        shell_file.write("needleall -asequence " + sample_ref+ " -bsequence "+i+ " -gapopen 10 -gapextend 0.5 -outfile sam/" + sample_name[:index_of_first_dash]+'_'+sample_name[index_of_first_dash+1:] +".sam -aformat sam &\n")
        
    else:
        pass
shell_file.write("wait\n")
shell_file.close()
