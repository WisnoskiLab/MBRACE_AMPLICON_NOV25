#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task 30
#SBATCH --job-name=Denoise

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.4

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs cutadapt/W_demux_seq_nov_2025_trim.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 280 \
  --o-table w_nov_2025_table.qza \
  --o-representative-sequences w_nov_2025_repseqs.qza \
  --o-denoising-stats w_nov_2025_denoising_stats.qza \
  --o-base-transition-stats w_nov_2025_base-transition-stats.qza