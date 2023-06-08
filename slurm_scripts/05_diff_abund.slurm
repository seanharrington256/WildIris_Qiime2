#!/bin/bash

#SBATCH --job-name diff_abund   # give the job a name
#SBATCH -A inbreh    # specify the account
#SBATCH -t 0-02:00		 # how much time? this is 1 hour
#SBATCH --nodes=1			# how many nodes?
#SBATCH --cpus-per-task=1	# 1 cores
#SBATCH --mem=10G			# 24 GB memory
#SBATCH --mail-type=ALL		# Send emails on start, fail, completion
#SBATCH --mail-user=EMAIL@EMAIL.com   # specify your email
#SBATCH -e errs_outs/err_diff_abund_%A.err		# name error files and include job ID (%A)
#SBATCH -o errs_outs/std_diff_abund_%A.out		# name standard out files and include job ID (%A)

# Load any necessary modules 
module load qiime2/2023.5

# set the working directory as appropriate
cd /project/wy_t3_2022/sharrin2/qiime_all_scripts/output 

# add psuedo-count:
qiime composition add-pseudocount \
  --i-table table.qza \
  --o-composition-table comp-table.qza


# ANCOM on vegetation
qiime composition ancom \
  --i-table comp-table.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization ancom-vegetation.qzv
  
# collapse features into sets of taxa, starting with the genus level
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-l6.qza

qiime composition add-pseudocount \
  --i-table table-l6.qza \
  --o-composition-table comp-table-l6.qza

qiime composition ancom \
  --i-table comp-table-l6.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization l6-ancom-vegetation.qzv
  
  
# Same for phyla:
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table table-l2.qza

qiime composition add-pseudocount \
  --i-table table-l2.qza \
  --o-composition-table comp-table-l2.qza

qiime composition ancom \
  --i-table comp-table-l2.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization l2-ancom-vegetation.qzv


  
  
