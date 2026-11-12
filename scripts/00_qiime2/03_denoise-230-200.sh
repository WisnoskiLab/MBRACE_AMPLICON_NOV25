#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task 60
#SBATCH --job-name=Denoise

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs cutadapt/W_demux_seq_nov_2025_trim.qza \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 200 \
  --o-table w_nov_2025_table_230_200.qza \
  --o-representative-sequences w_nov_2025_repseqs_230_200.qza \
  --o-denoising-stats w_nov_2025_denoising_stats_230_200.qza \
  --o-base-transition-stats w_nov_2025_base-transition-stats_230_200.qza