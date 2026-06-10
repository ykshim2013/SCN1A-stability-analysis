# Reproduction Guide

This guide documents the analysis pipeline. Per-variant input tables are **withheld until
publication**; the scripts expect them under a local `data/` directory. The source databases
(ClinVar, gnomAD v4.1, AlphaMissense, EVE) are public, so the calibration tables can be
reconstructed from them.

## Prerequisites

### Software
- **FoldX 5.1** — academic license from https://foldxsuite.crg.eu/
- **Rosetta** — academic license from https://www.rosettacommons.org/
- **Python 3.10+** with: pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl
  (see `requirements.txt`)

### Structure
The AlphaFold3-predicted Nav1.1 structure was generated with AlphaFold Server
(https://alphafoldserver.com/) from the full-length human SCN1A / Nav1.1 sequence
(UniProt **P35498**, 2,009 aa). The provided `structure/scn1a_wt_Repair.pdb` has already
been pre-processed with FoldX RepairPDB.

## Step 1 — FoldX stability calculations

```bash
python scripts/01_run_foldx_batch.py
```

Runs FoldX BuildModel with 5 independent iterations per variant. One directory per variant
is produced, each containing `Dif_*.fxout`, `Average_*.fxout`, `Raw_*.fxout`, and
`foldx_output.log`. Representative examples are provided in `raw_outputs/`.

## Step 2 — Parse FoldX results

```bash
python scripts/02_parse_foldx_results.py
```

Parses all `Dif_*.fxout` files and computes mean ΔΔG and SD across the 5 runs.

## Step 3 — Rosetta Cartesian ddG

```bash
cartesian_ddg.default.linuxgccrelease \
  -s structure/scn1a_wt_Repair.pdb \
  -ddg::mut_file mutations.mutfile \
  -ddg::iterations 3 \
  -score:weights ref2015_cart
```

Membrane-aware ΔΔG for the patch-clamp subset uses the `franklin2019` RosettaMP protocol.

## Step 4 — Conformational-state characterisation (optional)

```bash
python scripts/structural_state_diagnostic.py
```

Characterises the AF3 model's functional state (activation-gate minimum radius, DEKA
selectivity-filter geometry, per-domain S4 positions) for comparison with cryo-EM Nav
structures.

## Step 5 — Figures and downstream analysis

```bash
python scripts/06_generate_figures.py                 # stability / domain / topology figures
python scripts/07_population_predictor_analysis.py    # gnomAD-stratified + AlphaMissense/EVE comparison
```

These scripts also compute the reported discrimination statistics (ROC-AUC, Youden-optimal
threshold, domain-stratified enrichment).

## Expected behaviour (validation targets)

| Metric | Approximate value |
|--------|-------------------|
| FoldX ROC-AUC (P/LP vs benign) | ≈ 0.76 |
| Rosetta ROC-AUC | ≈ 0.70 |
| FoldX–Rosetta concordance (Spearman ρ) | ≈ 0.66 |
| Youden-optimal FoldX threshold | ΔΔG ≈ 1.60 kcal/mol |
| Specificity at ΔΔG > 2.0 kcal/mol | ≈ 0.98 |
| P-loop / selectivity-filter enrichment ratio | ≈ 0.40 |
| S4 voltage-sensor enrichment ratio | ≈ 2.32 |

## Computation environment

Analyses were run on macOS with FoldX 5.1, Rosetta (`ref2015_cart` / `franklin2019`),
Python 3.13, and the AlphaFold3 Server.
