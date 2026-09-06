#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=taxonomy

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.1

qiime feature-classifier classify-sklearn \
  --i-classifier SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
  --i-reads //mnt/scratch/wisnoskilab/dc2484/MBRACE/Sequences/November_2025/Qiime-2_analysis/W_qiime_work/w_nov_2025_repseqs.qza \
  --o-classification /mnt/scratch/wisnoskilab/dc2484/fce-lter/water_Sediment_Q2/Output_q2/w_nov_2025_taxonomy.qza