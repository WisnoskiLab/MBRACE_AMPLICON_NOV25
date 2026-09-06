#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task 20
#SBATCH --job-name=Import_Seq

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.4

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /mnt/scratch/wisnoskilab/dc2484/MBRACE/Sequences/November_2025/Qiime-2_analysis/S_qiime_work/W_manifest_nov_2025.tsv \
  --output-path W_demux_seq_nov_2025.qza \
  --input-format PairedEndFastqManifestPhred33V2