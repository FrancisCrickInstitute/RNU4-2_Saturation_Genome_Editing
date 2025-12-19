# RNU4-2_Saturation_Genome_Editing (SGE)

This is a repository containing the source code for analysing the RNU4-2 SGE data from De Jonghe _et al._ medrxiv (2025). Initial code from the pre-print released in April 2025 has now been supplemented with revision code. Revised scripts are contained in separate folders such that the original analysis on publicly available data can still be run prior to full publication of the article.


## Primary analysis
The primary analysis is performed on de-multiplexed fastq files (forward and reverse read files) and provides variant read count files for each replicate SGE experiment (e.g. RNU42brL41.txt and RNU42brL42.txt). The pipeline runs similarily to Buckley et _al._ Nature Genetics (2024). It comprises a main script RNU4-2_pipeline.sh which accesses a set of scripts to perform variant counting. It is recommended to generate a custom python environment with the following packages installed:
> Python/2.7.18-GCCcore-9.3.0

> GCC/5.4.0-2.26

> EMBOSS

The scripts and data files can be found in the primary_analysis directory, the original sequencing files corresponding to the data in the pre-print are available for download at ENA (project code).

## Secondary analysis
The secondary analysis directory takes variant read count files and performs filtering, normalisation and plotting. The main analysis is performed in a jupyter notebook (RNU4_2_SGE.ipynb). Some plots were generated in R and can be found in the RNU4_2_SGE.R script. All files to run the analysis can be found in the secondary_analysis directory.

A revised secondary analysis directory (secondary_analysis_v2) with code used in revision is now available. This contains revised scripts and additional scripts for new analyses added, such as the calibration of function scores to evidence codes. Input data tables (derived from new experiments done in revision of the manuscript) will be made publicly available upon publication.
