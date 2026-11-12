# MBRACE — November 2025 Estuarine 16S Amplicon Analysis

**Transport and Fate of Bacterial Communities in Mississippi Coastal Ecosystems**  
**Research group:** Wisnoski Lab, Mississippi State University  
**Sampling period:** November 2025  
**Locations:** Biloxi Bay and Pascagoula Bay, Mississippi, USA  
**Project background:** https://mbrace.usm.edu/transport-and-fate-of-bacterial-communities-in-mississippi-coastal-ecosystems

## Research scope

This repository organizes the November 2025 **16S rRNA gene amplicon** workflow for estuarine bacterial communities along freshwater–marine environmental gradients. The analyses investigate alpha diversity, beta diversity, taxonomic composition, spatial distance–decay, environmental drivers, abundance-defined subcommunities, niche width, and phylogenetic community assembly. It is an amplicon repository, not the master archive for MBRACE shotgun metagenomics or transcriptomics.

## Directory guide

```text
MBRACE_AMPLICON_NOV25/
├── README.md
├── .gitignore
├── MBRACE Analysis.Rproj
├── data/
│   ├── metadata/                # Sample metadata, environmental and station-distance tables
│   ├── processed/
│   │   ├── original/            # Earlier BIOM, FASTA, rooted/unrooted trees
│   │   └── trimmed_200_230/     # 200/230 bp variant used by current R scripts
│   └── derived/                 # Previously generated derived alpha table
├── scripts/
│   ├── 00_qiime2/               # HPC QIIME 2 pipeline scripts, ordered 01–05
│   ├── 01_alpha_beta/           # R code for diversity, salinity and NMDS
│   └── 02_spatial_community/    # R code for distance–decay, RDA, subcommunities, pNST
├── results/
│   ├── figures/
│   │   ├── alpha_diversity/
│   │   ├── beta_diversity/
│   │   ├── spatial_community/
│   │   └── legacy/             # Earlier plots retained for traceability
│   ├── tables/
│   │   ├── alpha_beta/
│   │   ├── spatial_community/
│   │   └── legacy/
│   └── models/                 # Saved model RDS objects
└── docs/
    ├── analysis_workflow.md
    ├── data_dictionary.md
    └── reorganization_notes.md
```

## Analysis sequence

### 1. QIIME 2 amplicon processing

Start with `scripts/00_qiime2/01_import.sh` followed by `02_cutadapt.sh`, `03_denoise-230-200.sh`, `04_align.sh`, and `05_taxa.sh`. These are SLURM scripts containing institutional paths, the QIIME environment name and specific sequence/parameter choices. **Check the lab HPC paths and parameters before running.** Raw FASTQ, QIIME artifacts and large intermediate files are not maintained as ordinary Git assets.

### 2. Alpha and beta diversity

`scripts/01_alpha_beta/Analysis-MBRACE.R` imports the processed 200/230 variant, representative sequences, phylogeny and metadata, constructs phyloseq/microeco objects, and produces richness, Shannon, salinity-group and NMDS analyses. Its input and output paths are repository-root relative. **Open `MBRACE Analysis.Rproj` or set the R working directory to the repository root.** The existing script uses a rarefaction depth of 40,000 reads; verify suitability and sample retention before treating results as final.

### 3. Spatial structure and assembly

`scripts/02_spatial_community/Distance-decay.R` covers abundance-defined subcommunities, geographic/environmental dissimilarity, RDA, niche width, and phylogenetic null-model/pNST analyses. Its known input paths were updated, but some model-cache and output paths remain relative to the execution directory and require review before rerunning. Pairwise comparisons are not independent replicates; confirm permutation procedures and statistical assumptions.

## File provenance and data variants

The repository previously held both `Analysis/Data files/` and `Analysis/Data files 200_230/`. Both processing variants have been retained separately, even where BIOM/taxon/tree compositions may differ. The original and 200/230 files must **not** be mixed without checking sample identifiers, representative sequences, ASV identifiers and tree tips. Byte-identical metadata/environmental duplicates were consolidated; non-identical files and distinct pNST model versions were retained. Refer to [data dictionary](docs/data_dictionary.md).

## Reproducibility

1. Clone the repository and open the project root in RStudio or VS Code.
2. Check `scripts/00_qiime2/` for HPC-specific environment and absolute paths.
3. Install the R libraries loaded by the chosen analysis script; dependency versions are not currently locked with `renv`.
4. Validate metadata, sample depths, ASV/tree correspondence and salinity group definitions.
5. Run scripts with the project root as working directory. Review relative outputs in the spatial script to avoid accidental writes outside `results/`.
6. Compare regenerated figures, tables and models with saved results before using them in a manuscript.

**Caution:** The reorganized repository preserves existing work; it does not indicate that scripts have been rerun or that statistical analyses have been independently validated.

## Project conventions

- Put code in `scripts/`, processed inputs in `data/`, figures in `results/figures/`, numerical outputs in `results/tables/` and cached models in `results/models/`.
- Make reproducible changes in feature branches and document major decisions in `docs/`.
- Do not commit OS files, R session state, credentials, new raw FASTQ files or scratch artifacts.
- Keep old versions clearly labeled until they can be scientifically evaluated.

## Use and attribution

Wisnoski Lab, Mississippi State University. For scientific reuse, authorship and citation, contact the project team. A public repository does not itself establish a reuse license.
