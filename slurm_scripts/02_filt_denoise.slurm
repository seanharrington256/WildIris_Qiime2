#!/bin/bash

#SBATCH --job-name filt_denoise   # give the job a name
#SBATCH -A inbreh    # specify the account
#SBATCH -t 0-02:00		 # how much time? this is 1 hour
#SBATCH --nodes=1			# how many nodes?
#SBATCH --cpus-per-task=1	# 1 cores
#SBATCH --mem=10G			# 10 GB memory
#SBATCH --mail-type=ALL		# Send emails on start, fail, completion
#SBATCH --mail-user=EMAIL@EMAIL.com   # specify your email
#SBATCH -e errs_outs/err_filt_denoise_%A.err		# name error files and include job ID (%A)
#SBATCH -o errs_outs/std_filt_denoise_%A.out		# name standard out files and include job ID (%A)

# Load any necessary modules 
module load qiime2/2023.5

# set the working directory as appropriate
cd /project/wy_t3_2022/sharrin2/qiime_all_scripts/output 

# Note that I changed the path to start in the output directory, since most of what we're working with is there now
# This changes how we specify the paths to input/output files below, mostly by just meaning we don't include
#  "output/" but also requiring us to put "../raw_data/" instead of just "raw_data"

# Filter samples with few reads
qiime tools export \
  --input-path demux-subsample.qzv \
  --output-path ./demux-subsample/

qiime demux filter-samples \
  --i-demux demux-subsample.qza \
  --m-metadata-file ./demux-subsample/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux demux.qza



# Denoise
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
  

# make some summaries:
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file ../raw_data/sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv  
