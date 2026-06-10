# Raw Computational Output Files

This directory contains unmodified output files from FoldX 5.1 and Rosetta Cartesian ddG for representative SCN1A variants. These files serve as verifiable evidence that stability calculations were genuinely performed using established structural biology software.

## FoldX Output Files

Each variant directory contains the following files produced by FoldX BuildModel:

### `Dif_scn1a_wt_Repair.fxout`
Energy difference between mutant and wild-type structures. Each row represents one independent run (5 total). Columns include total ddG and individual energy terms:
- Backbone Hbond, Sidechain Hbond, Van der Waals, Electrostatics
- Solvation Polar/Hydrophobic, Van der Waals clashes
- Entropy Sidechain/Mainchain, Torsional Clash, Helix Dipole
- Energy Ionisation, Partial Covalent Interactions

### `Average_scn1a_wt_Repair.fxout`
Averaged energy values across all 5 independent runs.

### `Raw_scn1a_wt_Repair.fxout`
Complete raw energy values for each individual run, showing both wild-type and mutant state energies.

### `foldx_output.log`
Execution log containing:
- FoldX 5.1 version header and license information
- Input parameter echo (numberOfRuns=5, pdb file, etc.)
- Per-run energy calculation details
- Completion status

### `individual_list.txt`
Mutation specification in FoldX format (e.g., `AA104D;` for variant A104D on chain A).

## Rosetta Output Files

### `mutations.ddg`
Rosetta Cartesian ddG output containing:
- COMPLEX energy scores for both WT and MUT states
- Three refinement rounds per state (Round 1, 2, 3)
- 18 individual Rosetta score function terms per round:
  `fa_atr, fa_rep, fa_sol, fa_intra_rep, lk_ball_wtd, fa_elec, hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, dslf_fa13, omega, fa_dun, p_aa_pp, yhh_planarity, ref, rama_prepro, cart_bonded`

### `mutations.mutfile`
Rosetta mutation specification format.

## Batch Run Logs

### `foldx_batch_865_run.log`
Complete log from the initial FoldX batch processing of 865 variants:
- Total time: 55.0 minutes
- 5 runs per variant, 8 parallel processes
- Per-variant status and ETA timestamps
- 865/865 successful, 0 failed

The 1,011-variant 5-run recalculation completed in 60.5 minutes (all 1,011 successful).

## Verification

### MD5 Checksums
`checksums.md5` contains MD5 hashes for all raw output files in this directory. Use these to verify file integrity:

```bash
# macOS
md5 -r <filename>

# Linux
md5sum <filename>
```

### Cross-reference with Manuscript
- A104D: FoldX ddG = 6.26 kcal/mol (check Dif file, average of 5 runs)
- R1648H: FoldX ddG = 5.93 kcal/mol (cited in manuscript Discussion)
- R101W: FoldX ddG = 14.32 kcal/mol (most destabilising variant; clipped at the y-axis maximum in the topology figure)
