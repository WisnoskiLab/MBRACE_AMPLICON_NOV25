# Reorganization notes — November 2026

- Split the former `Analysis/` directory into `data/`, `scripts/` and `results/`.
- Retained both processed ASV variants (`original/` and `trimmed_200_230/`).
- Consolidated byte-identical metadata/environmental copies; no distinct research dataset was intentionally discarded.
- Kept old plots in `results/figures/legacy` for traceability.
- Preserved the two distinct saved pNST model files as separate RDS objects.
- Removed only macOS `.DS_Store` files and adjusted known R input/output references.
- Some spatial script output paths and analysis settings still require review before execution.
- All changes should be checked before merging into `main`.
