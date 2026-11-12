#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task 20
#SBATCH --job-name=Import_Seq

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.4

mkdir -p cutadapt

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences W_demux_seq_nov_2025.qza \
  --p-front-f GTGYCAGCMGCCGCGGTA \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-error-rate 0.1 \
  --p-cores 20 \
  --o-trimmed-sequences cutadapt/W_demux_seq_nov_2025_trim.qza \
  --verbose

  qiime demux summarize \
  --i-data cutadapt/W_demux_seq_nov_2025_trim.qza \
  --o-visualization cutadapt/W_demux_seq_nov_2025_trim.qzv