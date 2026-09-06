#!/bin/bash  
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --job-name=merge 

qiime feature-classifier classify-sklearn \
  --i-classifier unite_ver2025-02-19_97_fungi-Q2-2026.4.qza \
  --i-reads representative_sequences.qza \
  --o-classification taxonomy.qza