#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task 20
#SBATCH --job-name=conver

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.4

qiime demux summarize --i-data W_demux_seq_nov_2025.qza --o-visualization W_demux-paired-end_nov_2025.qzv