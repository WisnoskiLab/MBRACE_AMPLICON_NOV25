#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task 20
#SBATCH --job-name=Import_Seq

conda activate qiime2-2026.4

FASTQ_DIR="/mnt/scratch/wisnoskilab/dc2484/MBRACE/Sequences/November_2025/Qiime-2_analysis/ITS"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$FASTQ_DIR" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ITS_demux_seq_nov_2025.qza
