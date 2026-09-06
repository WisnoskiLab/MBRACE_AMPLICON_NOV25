#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task 15
#SBATCH --job-name=Import_Seq

conda activate qiime2-2026.4

mkdir -p cutadapt

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences S_demux_seq_nov_2025.qza \
  --p-front-f GTGYCAGCMGCCGCGGTA \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-error-rate 0.1 \
  --p-cores 20 \
  --o-trimmed-sequences cutadapt/S_demux_seq_nov_2025_trim.qza \
  --verbose

  qiime demux summarize \
  --i-data cutadapt/S_demux_seq_nov_2025_trim.qza \
  --o-visualization cutadapt/S_demux_seq_nov_2025_trim.qzv