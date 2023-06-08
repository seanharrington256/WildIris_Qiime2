#!/bin/bash

#SBATCH --job-name imp_demux   # give the job a name
#SBATCH -A inbreh    # specify the account
#SBATCH -t 0-02:00		 # how much time? this is 1 hour
#SBATCH --nodes=1			# how many nodes?
#SBATCH --cpus-per-task=1	# 1 cores
#SBATCH --mem=10G			# 10 GB memory
#SBATCH --mail-type=ALL		# Send emails on start, fail, completion
#SBATCH --mail-user=EMAIL@EMAIL.com   # specify your email
#SBATCH -e errs_outs/err_import_%A.err		# name error files and include job ID (%A)
#SBATCH -o errs_outs/std_import_%A.out		# name standard out files and include job ID (%A)

# Load any necessary modules 
module load qiime2/2023.5

# set the working directory as appropriate
cd /project/wy_t3_2022/sharrin2/qiime_all_scripts

# Import the data into qiime format - read in from the raw_data directory and output into output directory
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path raw_data/emp-paired-end-sequences \
   --output-path output/emp-paired-end-sequences.qza


# Demultiplex
qiime demux emp-paired \
  --m-barcodes-file raw_data/sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs output/emp-paired-end-sequences.qza \
  --o-per-sample-sequences output/demux-full.qza \
  --o-error-correction-details output/demux-details.qza

# Subsample to make this tutorial faster
qiime demux subsample-paired \
  --i-sequences output/demux-full.qza \
  --p-fraction 0.3 \
  --o-subsampled-sequences output/demux-subsample.qza

# Summarize and make a visualization
qiime demux summarize \
  --i-data output/demux-subsample.qza \
  --o-visualization output/demux-subsample.qzv
