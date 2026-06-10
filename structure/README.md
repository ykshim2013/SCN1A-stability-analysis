# AlphaFold3 Structure Prediction

## Nav1.1 Wild-Type Structure

The full-length human Nav1.1 (SCN1A, UniProt P35498, 2,009 amino acids) structure was predicted using AlphaFold3 Server (https://alphafoldserver.com/).

## Quality Metrics

- **pTM:** 0.73
- **Ranking score:** 0.86
- **5 independent models** generated

## Per-Region Confidence (pLDDT)

| Region | Mean pLDDT | Classification |
|--------|-----------|----------------|
| S4 voltage sensors | 84.5 | Confident |
| VSD helices (S1-S3) | 87.1 | Confident |
| Pore helices (S5-S6) | 88.3 | Confident |
| P-loops | 89.6 | Confident |
| DI-DII linker | 37.0 | Very low |
| DII-DIII linker | 46.7 | Very low |
| C-terminus | 65.4 | Low |

## Files

- `scn1a_wt_Repair.pdb` — FoldX RepairPDB-optimized structure (used for all stability calculations)
- `alphafold3_confidence.json` — Per-residue pLDDT and pTM scores from model 0
- `alphafold3_job_request.json` — AlphaFold3 Server API request metadata

## Pre-processing

The AlphaFold3 structure was prepared using FoldX RepairPDB to optimize side-chain conformations prior to mutagenesis calculations. Global pre-relaxation was not performed.
