# Amplicon analysis workflow

1. **QIIME 2 processing:** `scripts/00_qiime2/01_import.sh` through `05_taxa.sh` (HPC-specific paths, modules and conda environment). Review parameters before execution.
2. **Diversity:** `scripts/01_alpha_beta/Analysis-MBRACE.R` (sample depth, rarefaction, group classification, alpha diversity and NMDS). Updated paths lead to the organized 200/230 data variant and results directories.
3. **Spatial and assembly:** `scripts/02_spatial_community/Distance-decay.R` (distance-decay, environment, RDA, niche width, pNST). Some legacy output filenames and cache paths remain inside the script; inspect before rerunning.

## Validation before publication

- Match BIOM ASV identifiers to representative sequences and phylogenetic tree tips.
- Confirm all sample metadata and environmental units and station codes.
- Check count filtering, singleton policies and rarefaction.
- Confirm transformations, permutational designs, random effects and bay-specific salinity support.
- Record software versions, random seeds, analysis settings and inputs.
- Rerun each script deliberately; files moved here are *previously generated* results, not newly validated outputs.
