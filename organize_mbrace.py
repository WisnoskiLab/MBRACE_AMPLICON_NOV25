#!/usr/bin/env python3
"""Non-destructive organization of WisnoskiLab/MBRACE_AMPLICON_NOV25.
Run at the repository root. Creates a git branch; preserves all scientifically distinct files.
"""
from pathlib import Path
import subprocess, shutil, filecmp

root = Path.cwd()
if not (root / '.git').exists() or not (root / 'Analysis').exists():
    raise SystemExit('Run this script from the cloned MBRACE_AMPLICON_NOV25 repository root (with Analysis/ and .git/ present).')

def git(*args):
    return subprocess.run(['git', *args], cwd=root, check=True, text=True, capture_output=True).stdout.strip()

if git('status', '--porcelain'):
    raise SystemExit('Git working tree has uncommitted changes. Commit or stash those changes first.')
branch = 'reorganize-amplicon-repository-2026'
if git('branch', '--list', branch):
    raise SystemExit(f'Branch {branch} already exists; inspect it before rerunning.')
git('switch', '-c', branch)

def move(src, dst):
    s, d = root / src, root / dst
    if not s.exists(): return
    d.parent.mkdir(parents=True, exist_ok=True)
    if d.exists():
        if s.is_file() and d.is_file() and filecmp.cmp(s,d,shallow=False):
            s.unlink(); print('Deduplicated',src); return
        raise SystemExit(f'Conflict: {src} and {dst} are distinct. Resolve manually.')
    shutil.move(str(s), str(d)); print(f'{src} -> {dst}')

# Move shared metadata first, then consolidate duplicates using byte comparison.
for name in ['metadata.txt','nov_2025_env_data copy.xlsx','nov_2025_env_data.txt','nov_2025_env_data.xlsx','station_transect_distances.csv']:
    dest = {'metadata.txt':'metadata_200_230.txt','nov_2025_env_data copy.xlsx':'nov_2025_env_data_copy.xlsx'}.get(name,name)
    move('Analysis/Data files 200_230/'+name, 'data/metadata/'+dest)
for name in ['metadata.txt','nov_2025_env_data copy.xlsx','nov_2025_env_data.txt','nov_2025_env_data.xlsx','station_transect_distances.csv']:
    dest = {'metadata.txt':'metadata_original.txt','nov_2025_env_data copy.xlsx':'nov_2025_env_data_copy.xlsx'}.get(name,name)
    move('Analysis/Data files/'+name, 'data/metadata/'+dest)
for name in ['metadata_nov_2025.xlsx','mb_analyst_metadata.csv']:
    move('Analysis/Data files/'+name, 'data/metadata/'+name)
move('Analysis/Data files/Alpha-edited.csv','data/derived/Alpha-edited.csv')
for name in ['w_nov_2025-rooted_tree.nwk','w_nov_2025-unrooted_tree.nwk','w_nov_2025_repseqs.fasta','w_nov_2025_table.biom']:
    move('Analysis/Data files/'+name,'data/processed/original/'+name)
for name in ['dna-sequences_230_200.fasta','w_nov_2025_200-230_feture.biom','w_nov_2025_200_230_rooted-tree.nwk','w_nov_2025_200_230_unrooted-tree.nwk']:
    move('Analysis/Data files 200_230/'+name,'data/processed/trimmed_200_230/'+name)
for name in ['01_import.sh','02_cutadapt.sh','03_denoise-230-200.sh','04_align.sh','05_taxa.sh']:
    move('Analysis/Qiime Scripts/'+name,'scripts/00_qiime2/'+name)
move('Analysis/Analysis-MBRACE.R','scripts/01_alpha_beta/Analysis-MBRACE.R')
move('Analysis/Distance-decay.R','scripts/02_spatial_community/Distance-decay.R')
move('Analysis/pNST_models_corrected.rds','results/models/pNST_models_corrected_previous.rds')
for p in sorted((root/'Analysis/Distance-Decay results').glob('*')):
    if not p.is_file(): continue
    dest='results/models' if p.suffix=='.rds' else ('results/tables/spatial_community' if p.suffix=='.csv' else 'results/figures/spatial_community')
    move(str(p.relative_to(root)),dest+'/'+p.name)
for p in sorted((root/'Analysis/Output').rglob('*')):
    if not p.is_file() or p.name=='.DS_Store': continue
    part=p.parent.name
    dest={'Alpha':'results/figures/alpha_diversity','Beta':'results/figures/beta_diversity','Tables':'results/tables/alpha_beta'}.get(part)
    if dest: move(str(p.relative_to(root)),dest+'/'+p.name)
for p in sorted((root/'Analysis/Output Results').glob('*')):
    if p.is_file() and p.name!='.DS_Store':
        move(str(p.relative_to(root)),('results/figures/legacy' if p.suffix=='.pdf' else 'results/tables/legacy')+'/'+p.name)

# Adjust the *known* file inputs and outputs. No numerical analysis is changed.
a = root/'scripts/01_alpha_beta/Analysis-MBRACE.R'
if a.exists():
    s=a.read_text().replace('Data files 200_230/','data/processed/trimmed_200_230/').replace('Output/Alpha','results/figures/alpha_diversity').replace('Output/Beta','results/figures/beta_diversity').replace('Output/Tables','results/tables/alpha_beta')
    a.write_text('# Run from the repository root.\n'+s)
b = root/'scripts/02_spatial_community/Distance-decay.R'
if b.exists():
    s=b.read_text().replace('Data files 200_230/','data/processed/trimmed_200_230/').replace('Data files/nov_2025_env_data copy.xlsx','data/metadata/nov_2025_env_data_copy.xlsx').replace('Data files/station_transect_distances.csv','data/metadata/station_transect_distances.csv')
    b.write_text('# Run from repository root; review remaining output paths before rerunning.\n'+s)

for p in [root/'.DS_Store',*root.glob('Analysis/**/.DS_Store')]:
    if p.is_file(): p.unlink()
for p in sorted((root/'Analysis').rglob('*'), key=lambda x:len(x.parts), reverse=True):
    if p.is_dir() and not any(p.iterdir()): p.rmdir()
if (root/'Analysis').is_dir() and not any((root/'Analysis').iterdir()): (root/'Analysis').rmdir()

(root/'.gitignore').write_text('''.Rproj.user/\n.Rhistory\n.RData\n.Ruserdata\n.DS_Store\n**/.DS_Store\n.Renviron\n.env\n.env.*\n*.fastq\n*.fastq.gz\n*.fq\n*.fq.gz\n*.qza\n*.qzv\n*.tmp\n''')
readme = '''# MBRACE — November 2025 Estuarine 16S Amplicon Analysis

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
'''
(root/'README.md').write_text(readme)
(root/'docs').mkdir(exist_ok=True)
(root/'docs'/'analysis_workflow.md').write_text('''# Amplicon analysis workflow

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
''')
(root/'docs'/'data_dictionary.md').write_text('''# Data dictionary: file-level inventory

| Path | Contents |
|---|---|
| `data/metadata/metadata_200_230.txt` | Metadata packaged with the 200/230 export |
| `data/metadata/metadata_original.txt` | Same-content metadata packaged with the original export |
| `data/metadata/metadata_nov_2025.xlsx` | Additional sample metadata workbook |
| `data/metadata/mb_analyst_metadata.csv` | Metadata adapted for another analysis workflow |
| `data/metadata/nov_2025_env_data.xlsx` | Environmental measurements workbook |
| `data/metadata/nov_2025_env_data_copy.xlsx` | Alternate environmental workbook, retained separately |
| `data/metadata/nov_2025_env_data.txt` | Environmental text export |
| `data/metadata/station_transect_distances.csv` | Transect/station distance data |
| `data/processed/original/` | Original 16S BIOM table, representative sequences and phylogenies |
| `data/processed/trimmed_200_230/` | 200/230 processed 16S BIOM, sequences and phylogenies |
| `data/derived/Alpha-edited.csv` | Existing derived alpha-diversity table |

**TODO:** Validate the headers and units of each metadata column, the exact sample identifiers and the QIIME 2/reference database versions against lab records before formal data release.
''')
(root/'docs'/'reorganization_notes.md').write_text('''# Reorganization notes — November 2026

- Split the former `Analysis/` directory into `data/`, `scripts/` and `results/`.
- Retained both processed ASV variants (`original/` and `trimmed_200_230/`).
- Consolidated byte-identical metadata/environmental copies; no distinct research dataset was intentionally discarded.
- Kept old plots in `results/figures/legacy` for traceability.
- Preserved the two distinct saved pNST model files as separate RDS objects.
- Removed only macOS `.DS_Store` files and adjusted known R input/output references.
- Some spatial script output paths and analysis settings still require review before execution.
- All changes should be checked before merging into `main`.
''')
print('\nFinished on branch:', branch)
print('Review changes using: git status --short && git diff --stat && git diff -- README.md')
print('Then: git add -A && git commit -m "Reorganize MBRACE amplicon analysis repository" && git push -u origin '+branch)
print('Create a pull request from the branch into main on GitHub.')
