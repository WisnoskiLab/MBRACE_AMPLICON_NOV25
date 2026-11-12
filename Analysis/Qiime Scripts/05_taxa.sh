#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=60
#SBATCH --job-name=taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza \
  --i-reads w_nov_2025_repseqs_230_200.qza \
  --o-classification w_nov_2025_taxonomy230_200.qza