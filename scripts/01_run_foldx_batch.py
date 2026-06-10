#!/usr/bin/env python3
"""
FoldX Batch Processing - All 1,011 SCN1A variants with 5 runs each.
Recalculation for peer review revision (reviewer requested 5 iterations).
"""

import os
import subprocess
import shutil
import time
import re
from multiprocessing import Pool, cpu_count
from pathlib import Path
import pandas as pd

# Configuration
FOLDX = "/Users/ykshim2025/Desktop/Code2025/FoldX5.1/foldx5_Mac/foldx_20261231"
WORKDIR = Path("/Users/ykshim2025/Desktop/Code2025/SCN1A/foldx_analysis")
PDB_FILE = "scn1a_wt_Repair.pdb"
OUTPUT_DIR = WORKDIR / "foldx_results_1011_5runs"
VARIANTS_FILE = WORKDIR / "scn1a_variants_1011.txt"
NUM_CORES = min(cpu_count(), 8)
NUM_RUNS = 5


def run_foldx_variant(variant):
    """Run FoldX BuildModel for a single variant with 5 runs"""
    variant_dir = OUTPUT_DIR / variant

    try:
        # Skip if already successfully completed
        dif_file = list(variant_dir.glob("Dif_*.fxout")) if variant_dir.exists() else []
        if dif_file:
            return variant, "SKIPPED", None

        variant_dir.mkdir(parents=True, exist_ok=True)

        # Copy PDB
        pdb_path = variant_dir / PDB_FILE
        if not pdb_path.exists():
            shutil.copy(WORKDIR / PDB_FILE, pdb_path)

        # Parse variant name (e.g., A104D -> AA104D;)
        match = re.match(r'^([A-Z])(\d+)([A-Z])$', variant)
        if not match:
            return variant, "ERROR", f"Could not parse variant: {variant}"

        wt_aa, pos, mut_aa = match.groups()
        foldx_mutation = f"{wt_aa}A{pos}{mut_aa};"

        individual_list = variant_dir / "individual_list.txt"
        with open(individual_list, 'w') as f:
            f.write(foldx_mutation + "\n")

        cmd = [
            FOLDX,
            "--command=BuildModel",
            f"--pdb={PDB_FILE}",
            "--mutant-file=individual_list.txt",
            f"--numberOfRuns={NUM_RUNS}",
            "--out-pdb=false"
        ]

        with open(variant_dir / "foldx_output.log", 'w') as logfile:
            result = subprocess.run(
                cmd,
                cwd=str(variant_dir),
                stdout=logfile,
                stderr=subprocess.STDOUT,
                timeout=3600
            )

        if result.returncode == 0:
            return variant, "SUCCESS", None
        else:
            return variant, "FAILED", f"Return code: {result.returncode}"

    except subprocess.TimeoutExpired:
        return variant, "TIMEOUT", "Exceeded 1 hour timeout"
    except Exception as e:
        return variant, "ERROR", str(e)


def parse_results():
    """Parse all FoldX Dif_ output files into a summary TSV"""
    records = []
    for variant_dir in sorted(OUTPUT_DIR.iterdir()):
        if not variant_dir.is_dir():
            continue
        variant = variant_dir.name
        dif_files = list(variant_dir.glob("Dif_*.fxout"))
        if not dif_files:
            continue

        dif_file = dif_files[0]
        ddg_values = []
        with open(dif_file) as f:
            for line in f:
                if line.startswith("scn1a_wt_Repair") or line.startswith(PDB_FILE.replace(".pdb", "")):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        try:
                            ddg_values.append(float(parts[1]))
                        except ValueError:
                            pass

        if ddg_values:
            import numpy as np
            records.append({
                'Variant': variant,
                'FoldX_ddG_5runs': np.mean(ddg_values),
                'FoldX_ddG_std_5runs': np.std(ddg_values, ddof=1) if len(ddg_values) > 1 else 0.0,
                'n_runs': len(ddg_values)
            })

    if records:
        df = pd.DataFrame(records)
        out_path = OUTPUT_DIR / "foldx_ddg_1011_5runs.tsv"
        df.to_csv(out_path, sep='\t', index=False)
        print(f"\nParsed {len(df)} variants -> {out_path}")
        print(f"  n_runs distribution: {df['n_runs'].value_counts().to_dict()}")
        print(f"  Median SD: {df['FoldX_ddG_std_5runs'].median():.4f}")
        print(f"  Mean ddG: {df['FoldX_ddG_5runs'].mean():.3f}")
        return df
    return None


def main():
    print("=" * 70)
    print("FoldX 5-Run Recalculation - 1,011 SCN1A Variants")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(VARIANTS_FILE) as f:
        variants = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    print(f"Total variants: {len(variants)}")
    print(f"Runs per variant: {NUM_RUNS}")
    print(f"Parallel processes: {NUM_CORES}")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 70)

    start_time = time.time()
    success = 0
    failed = 0
    skipped = 0
    failed_list = []

    with Pool(NUM_CORES) as pool:
        for i, result in enumerate(pool.imap_unordered(run_foldx_variant, variants)):
            variant, status, error = result
            elapsed = time.time() - start_time
            done = i + 1
            eta = (elapsed / done) * (len(variants) - done)

            if status == "SUCCESS":
                success += 1
                print(f"[{done}/{len(variants)}] {variant}: SUCCESS (ETA: {eta/60:.1f} min)")
            elif status == "SKIPPED":
                skipped += 1
                if done % 100 == 0:
                    print(f"[{done}/{len(variants)}] ... {skipped} skipped so far (ETA: {eta/60:.1f} min)")
            else:
                failed += 1
                failed_list.append((variant, status, error))
                print(f"[{done}/{len(variants)}] {variant}: {status} - {error}")

    total_time = time.time() - start_time

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total: {len(variants)} | Success: {success} | Skipped: {skipped} | Failed: {failed}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print("=" * 70)

    if failed_list:
        print("\nFailed variants:")
        for v, s, e in failed_list:
            print(f"  {v}: {s} - {e}")

    # Parse results
    print("\nParsing FoldX output files...")
    parse_results()


if __name__ == "__main__":
    main()
