#!/bin/bash

#SBATCH --job-name phylo_div   # give the job a name
#SBATCH -A inbreh    # specify the account
#SBATCH -t 0-02:00		 # how much time? this is 1 hour
#SBATCH --nodes=1			# how many nodes?
#SBATCH --cpus-per-task=1	# 1 cores
#SBATCH --mem=10G			# 10 GB memory
#SBATCH --mail-type=ALL		# Send emails on start, fail, completion
#SBATCH --mail-user=EMAIL@EMAIL.com   # specify your email
#SBATCH -e errs_outs/err_phylo_div_%A.err		# name error files and include job ID (%A)
#SBATCH -o errs_outs/std_phylo_div_%A.out		# name standard out files and include job ID (%A)

# Load any necessary modules 
module load qiime2/2023.5

# set the working directory as appropriate
cd /project/wy_t3_2022/sharrin2/qiime_all_scripts/output 

# Make a phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
# Calculate diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 689 \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --output-dir core-metrics-results  
  
  
  
# categorical tests using the evenness and Faith Phylogenetic Diversity metrics
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
  

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv  

# Test if beta diversity varies with vegetation
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --m-metadata-column vegetation \
  --o-visualization core-metrics-results/unweighted-unifrac-vegetation-significance.qzv \
  --p-pairwise



# Alpha diversity correlations
qiime diversity alpha-correlation \
	--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	--m-metadata-file ../raw_data/sample-metadata.tsv \
	--p-method spearman \
	--o-visualization core-metrics-results/faith-pd-correlation.qzv




# For beta diversity correlations, first calculate distance matrix
qiime metadata distance-matrix \
	--m-metadata-file ../raw_data/sample-metadata.tsv \
	--m-metadata-column elevation \
	--o-distance-matrix core-metrics-results/elevation-dist-mat.qza

# Mantel test to test for an association between this elevation distance matrix and one of our metrics of beta diversity
qiime diversity mantel \
	--i-dm1 core-metrics-results/elevation-dist-mat.qza \
	--i-dm2 core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--p-method spearman \
	--p-label1 elevation \
	--p-label2 unweighted_unifrac \
	--p-intersect-ids \
	--o-visualization core-metrics-results/unweight_unifrac_elevation_mantel.qzv


# Alpha rarefaction:
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1175 \
  --m-metadata-file ../raw_data/sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv




