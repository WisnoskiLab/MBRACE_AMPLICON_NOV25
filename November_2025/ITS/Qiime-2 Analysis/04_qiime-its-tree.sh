#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --job-name=ITS_tree
#SBATCH --output=ITS_tree_%j.out

set -euo pipefail

echo "========================================"
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "========================================"

# Activate your QIIME2 environment here
# Example:
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate qiime2-amplicon-XXXX

echo "QIIME location:"
which qiime

echo "QIIME version:"
qiime info

echo "Checking input file:"
ls -lh representative_sequences.qza

echo "Starting phylogenetic tree construction..."
date

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative_sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza \
  --verbose

echo "Tree construction finished."
date

echo "Output files:"
ls -lh \
  aligned-rep-seqs.qza \
  masked-aligned-rep-seqs.qza \
  unrooted-tree.qza \
  rooted-tree.qza

echo "Job completed: $(date)"