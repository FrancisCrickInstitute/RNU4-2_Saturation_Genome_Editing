#RNU4_2_b.sh

#This updated analysis works off of raw reads, starting with bcl2fastq ; also moves over files from asf with fastqc data, too - redundant for now
#where is sequencing data?

#Change PATH with running PATH
out_dir="RNU4_2_b"

adapter1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" #this is the nextera reverse adapter, as seen in R1
adapter2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"  #this is the nextera forward adapter, as seen in R2
seq_type="N" #or Miseq, change to N if nextseq
ref_folder="RNU4_2_b/fasta"
#must not include spaces:
cigar_comparisons="RNU42rL42_D5+RNU42rL42_D14"
amplicon_list="RNU42" #splits on commas
#6 entries per experimental grouping must be specified, order must be:  experiment, pre, post, lib, neg, rna, as in here:
exp_groupings="RNU42rL41+RNU42rL41_D5+RNU42rL41_D14+RNU42_lib+RNU42_neg+RNU42rL41_D14,RNU42rL42+RNU42rL42_D5+RNU42rL42_D14+RNU42_lib+RNU42_neg+RNU42rL42_D14" #splits on commas, then +
#A file listing the expected pam edits, the sge region in the amplicon, and which gRNA was used, etc.
editing_info_file="RNU42b_editing_info.txt"

#for determining whether or not to include a read in output of sam_to_edits
#reads_threshold=".000002"
reads_threshold=".000002"
alignment_score_threshold="300"

#Is gene on "+" or "-" strand of genome (e.g. "-" for BRCA1, "+" for VHL)
orientation="-"
#Consensus coding sequence for gene, should match what's used in ClinVar
CCDS="CCDS33356.1"
#A tab-delimited text file listing the amplicons and their genomic coordinates (must match CADD human genome version, e.g. hg19 or hg38)
amplicon_coords="RNU42_amp_coords.txt"
#directory with cadd files for each amplicon, named as {amplicon}.cadd 
min_indel_freq=".00002"

ml load Python/2.7.18-GCCcore-9.3.0

mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/PATH/run_seqprep_210729.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt
mkdir no_Ns
python ~/PATH/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns

ml load Python/2.7.18-GCCcore-9.3.0
module load GCC/5.4.0-2.26
module load EMBOSS


mkdir sam
python ~/PATH/run_needle_to_sam_RNU42_pipeline.py $ref_folder
echo "Running needleall."

ml purge
sh run_needle_to_sam.sh #input dashes, output is underscored
wait #lunch
echo "Finished running needleall."
cd sam
mkdir cigar_counts

python ~/PATH/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -15 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt
mkdir variant_counts
#writes minimally thresholded output (set in hardcode (0.000001 in any) to variant_counts folder
python 210317_sam_to_edits_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold