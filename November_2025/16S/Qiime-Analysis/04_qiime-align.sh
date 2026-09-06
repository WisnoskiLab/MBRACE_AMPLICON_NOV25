#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=200G
#SBATCH --cpus-per-task 20
#SBATCH --job-name=align

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate

conda activate qiime2-2026.4

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences w_nov_2025_repseqs.qza \
  --o-alignment w_nov_2025_aligned-rep-seqs.qza \
  --o-masked-alignment w_nov_2025_masked-aligned-rep-seqs.qza \
  --o-tree w_nov_2025_unrooted-tree.qza \
  --o-rooted-tree w_nov_2025_rooted-tree.qza