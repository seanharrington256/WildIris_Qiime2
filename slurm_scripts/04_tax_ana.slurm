#!/bin/bash

#SBATCH --job-name tax_ana   # give the job a name
#SBATCH -A inbreh    # specify the account
#SBATCH -t 0-02:00		 # how much time? this is 1 hour
#SBATCH --nodes=1			# how many nodes?
#SBATCH --cpus-per-task=1	# 1 cores
#SBATCH --mem=24G			# 24 GB memory
#SBATCH --mail-type=ALL		# Send emails on start, fail, completion
#SBATCH --mail-user=EMAIL@EMAIL.com   # specify your email
#SBATCH -e errs_outs/err_tax_ana_%A.err		# name error files and include job ID (%A)
#SBATCH -o errs_outs/std_tax_ana_%A.out		# name standard out files and include job ID (%A)

# Load any necessary modules 
module load qiime2/2023.5

# set the working directory as appropriate
cd /project/wy_t3_2022/sharrin2/qiime_all_scripts/output 



# extract reads
qiime feature-classifier extract-reads \
  --i-sequences $QIIME_DATA_RESOURCES/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy $QIIME_DATA_RESOURCES/silva-138-99-tax.qza \
  --o-classifier classifier.qza


# Classify our data:
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


# tabulate the output
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Make a barplot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

