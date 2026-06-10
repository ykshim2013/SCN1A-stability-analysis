# SCN1A Protein Stability Analysis

Computational protein-stability methodology for SCN1A (Nav1.1) missense variants, combining FoldX and Rosetta Cartesian ΔΔG on an AlphaFold3-predicted Nav1.1 structure to distinguish a proteostasis-impaired, pore-region variant subset from a complementary stability-neutral population.

## Publication

**Title:** Computational Protein Stability Analysis of SCN1A Missense Variants Identifies a Proteostasis-Impaired Pore-Region Subset and a Complementary Stability-Neutral Population

**Journal:** *Epilepsia Open* (under review)

> **Scope of this repository.** This repository documents the **methodology and analytical concept** of the study. The manuscript, full per-variant datasets, computed result tables, and publication figures are **withheld until the paper is accepted**. The analysis code, the AlphaFold3 input structure, representative raw tool outputs, and the reproduction guide are provided so the pipeline can be inspected and re-run against the public source databases.

## Concept

Missense variants in *SCN1A* can cause epilepsy through at least two structurally distinct routes: (1) **protein destabilisation** that promotes misfolding, ER retention, and proteostatic degradation (reducing functional channel density), and (2) **gating dysfunction** in channels that remain structurally stable. This study uses physics-based ΔΔG prediction on a full-length Nav1.1 model to estimate per-variant destabilisation, then stratifies the signal by functional domain to identify which variant subsets are most consistent with each mechanism.

## Key Findings

- Predicted destabilisation is significantly associated with pathogenicity (FoldX ROC-AUC ≈ 0.76; Rosetta Cartesian ΔΔG concordant, Spearman ρ = 0.66).
- **P-loop / selectivity-filter** pathogenic variants are strongly depleted of stability-neutral cases (0.40-fold; Bonferroni-adjusted *p* = 4.3 × 10⁻⁷) — consistent with a proteostasis-impaired interpretation.
- **S4 voltage-sensor** pathogenic variants are enriched for stability-neutral cases (2.32-fold; adjusted *p* = 0.013) — a hypothesis-generating, gating-dysfunction–leaning subset.
- Strong destabilisation (ΔΔG > 2.0 kcal/mol) gives high specificity (≈ 0.98) as supporting evidence for pathogenicity; the Youden-optimal discrimination threshold is ΔΔG = 1.60 kcal/mol.
- Stability prediction is **complementary** to sequence/conservation predictors (AlphaMissense, EVE), which achieve higher overall discrimination but do not resolve mechanism.

## Computational Methods

| Tool | Version / protocol | Iterations | Purpose |
|------|--------------------|------------|---------|
| AlphaFold3 | AlphaFold Server | 5 models | Full-length Nav1.1 structure prediction |
| FoldX | 5.1 — RepairPDB + BuildModel | 5 per variant | Stability prediction (ΔΔG, kcal/mol) |
| Rosetta | Cartesian ddG (`ref2015_cart`) | 3 per variant | Independent ΔΔG validation (REU) |
| RosettaMP | `franklin2019` | 3 per variant | Membrane-aware stability (patch-clamp subset) |

Variant references: ClinVar (P/LP, B/LB, VUS) and gnomAD v4.1 (allele-frequency–stratified population comparator). External predictors: AlphaMissense and EVE.

## Repository Structure

```
scripts/            Analysis pipeline (FoldX run/parse, figures, population & predictor
                    analysis, conformational-state diagnostic)
structure/          AlphaFold3-predicted Nav1.1 structure + per-residue confidence
raw_outputs/        Representative unmodified FoldX/Rosetta outputs (3 example variants)
reproducibility/    Environment specification and step-by-step reproduction guide
LICENSE             CC-BY 4.0
```

> Input datasets (`data/`), computed results (`results/`), supplementary tables, and
> rendered figures are intentionally **not** included pre-acceptance. The scripts expect
> input tables under a local `data/` directory; see `reproducibility/REPRODUCE.md`.

## Scripts

| Script | Role |
|--------|------|
| `scripts/01_run_foldx_batch.py` | FoldX BuildModel (5-run) batch over all variants |
| `scripts/02_parse_foldx_results.py` | Parse `Dif_*.fxout`, aggregate mean/SD ΔΔG |
| `scripts/06_generate_figures.py` | Generate stability / domain / topology figures |
| `scripts/07_population_predictor_analysis.py` | gnomAD-stratified and AlphaMissense/EVE comparison analysis |
| `scripts/structural_state_diagnostic.py` | Conformational-state characterisation of the AF3 model (gate radius, DEKA geometry, S4 positions) |

## Raw Computational Evidence

`raw_outputs/` contains unmodified FoldX and Rosetta outputs for representative variants
(A104D, R101W, R1648H) — `.fxout`/`.ddg` files with full energy-term breakdowns, execution
logs, and MD5 checksums — as verifiable evidence that the calculations were genuinely run
with established structural-biology software. See [`raw_outputs/README.md`](raw_outputs/README.md).

## Reproduction

See [`reproducibility/REPRODUCE.md`](reproducibility/REPRODUCE.md) for prerequisites and the
step-by-step pipeline.

## Data Availability

- **ClinVar:** https://www.ncbi.nlm.nih.gov/clinvar/
- **gnomAD v4.1:** https://gnomad.broadinstitute.org/
- **AlphaFold3 Server:** https://alphafoldserver.com/
- **AlphaMissense / EVE:** respective public releases

The complete per-variant dataset and publication figures will be released here upon acceptance.

## License

Licensed under CC-BY 4.0 (see `LICENSE`).

## Citation

[To be updated upon publication]
