#!/usr/bin/env python3
"""
Parse FoldX 5-Run Results - Publication Standard Analysis
Calculates mean and median ddG from 5 independent runs
"""

from pathlib import Path
import statistics

WORKDIR = Path("/Users/ykshim2025/Desktop/Code2025/SCN1A/foldx_analysis")
RESULTS_DIR = WORKDIR / "foldx_results_5runs"
OUTPUT_FILE = RESULTS_DIR / "foldx_ddg_summary_5runs.tsv"

# Classification mapping based on variant lists
PATHOGENIC_ORIGINAL = set(open(WORKDIR / "scn1a_variants.txt").read().split()) - {'#'}
PATHOGENIC_NEW = set(open(WORKDIR / "scn1a_variants_new.txt").read().split()) - {'#'}
PATHOGENIC_ADDITIONAL = set(open(WORKDIR / "scn1a_variants_additional.txt").read().split()) - {'#'}
BENIGN = set(open(WORKDIR / "scn1a_variants_benign.txt").read().split()) - {'#'}

def get_classification(variant):
    """Determine variant classification"""
    # Clean variant name
    v = variant.strip()
    if v in BENIGN or v.startswith('#'):
        return "Benign"
    return "Pathogenic"

def get_source(variant):
    """Determine variant source"""
    v = variant.strip()
    if v in BENIGN:
        return "ClinVar-Benign"
    if v in PATHOGENIC_ADDITIONAL:
        return "Additional-GOF"
    if v in PATHOGENIC_NEW:
        return "Expanded-2024"
    return "Original"

def parse_dif_file(filepath):
    """Parse FoldX Dif_*.fxout file to get all ddG values from 5 runs"""
    ddg_values = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('scn1a_wt_Repair'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    try:
                        ddg = float(parts[1])
                        ddg_values.append(ddg)
                    except ValueError:
                        continue
    return ddg_values

def classify_effect(ddg):
    """Classify stability effect based on ddG"""
    if ddg > 2:
        return "Highly Destabilizing"
    elif ddg > 0.5:
        return "Destabilizing"
    elif ddg > -0.5:
        return "Neutral"
    elif ddg > -2:
        return "Stabilizing"
    else:
        return "Highly Stabilizing"

def main():
    results = []

    for variant_dir in sorted(RESULTS_DIR.iterdir()):
        if not variant_dir.is_dir():
            continue
        variant = variant_dir.name
        dif_file = variant_dir / "Dif_scn1a_wt_Repair.fxout"

        if dif_file.exists():
            ddg_values = parse_dif_file(dif_file)
            if ddg_values:
                mean_ddg = statistics.mean(ddg_values)
                median_ddg = statistics.median(ddg_values)
                std_ddg = statistics.stdev(ddg_values) if len(ddg_values) > 1 else 0

                results.append({
                    'variant': variant,
                    'ddG_mean': mean_ddg,
                    'ddG_median': median_ddg,
                    'ddG_std': std_ddg,
                    'n_runs': len(ddg_values),
                    'effect': classify_effect(mean_ddg),
                    'classification': get_classification(variant),
                    'source': get_source(variant)
                })
            else:
                results.append({
                    'variant': variant,
                    'ddG_mean': None,
                    'ddG_median': None,
                    'ddG_std': None,
                    'n_runs': 0,
                    'effect': 'Parse Error',
                    'classification': get_classification(variant),
                    'source': get_source(variant)
                })
        else:
            results.append({
                'variant': variant,
                'ddG_mean': None,
                'ddG_median': None,
                'ddG_std': None,
                'n_runs': 0,
                'effect': 'Failed',
                'classification': get_classification(variant),
                'source': get_source(variant)
            })

    # Write detailed results
    with open(OUTPUT_FILE, 'w') as f:
        f.write("Variant\tddG_mean\tddG_median\tddG_std\tN_Runs\tEffect\tClassification\tSource\n")
        for r in results:
            if r['ddG_mean'] is not None:
                f.write(f"{r['variant']}\t{r['ddG_mean']:.3f}\t{r['ddG_median']:.3f}\t{r['ddG_std']:.3f}\t{r['n_runs']}\t{r['effect']}\t{r['classification']}\t{r['source']}\n")
            else:
                f.write(f"{r['variant']}\tNA\tNA\tNA\t{r['n_runs']}\t{r['effect']}\t{r['classification']}\t{r['source']}\n")

    # Statistical analysis
    successful = [r for r in results if r['ddG_mean'] is not None]
    pathogenic = [r for r in successful if r['classification'] == 'Pathogenic']
    benign = [r for r in successful if r['classification'] == 'Benign']

    path_ddg = [r['ddG_mean'] for r in pathogenic]
    benign_ddg = [r['ddG_mean'] for r in benign]

    print("=" * 70)
    print("FoldX 5-Run Results Summary (Publication Standard)")
    print("=" * 70)
    print(f"Total analyzed: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(results) - len(successful)}")
    print(f"Runs per variant: 5")

    print("\n" + "-" * 70)
    print("PATHOGENIC VARIANTS")
    print("-" * 70)
    print(f"N: {len(pathogenic)}")
    if path_ddg:
        print(f"Mean ddG: {statistics.mean(path_ddg):.2f} +/- {statistics.stdev(path_ddg):.2f} kcal/mol")
        print(f"Median ddG: {statistics.median(path_ddg):.2f} kcal/mol")
        print(f"Range: [{min(path_ddg):.2f}, {max(path_ddg):.2f}] kcal/mol")

        effects = {}
        for r in pathogenic:
            effects[r['effect']] = effects.get(r['effect'], 0) + 1
        print("Effect distribution:", effects)

    print("\n" + "-" * 70)
    print("BENIGN VARIANTS (Control)")
    print("-" * 70)
    print(f"N: {len(benign)}")
    if benign_ddg:
        print(f"Mean ddG: {statistics.mean(benign_ddg):.2f} +/- {statistics.stdev(benign_ddg):.2f} kcal/mol")
        print(f"Median ddG: {statistics.median(benign_ddg):.2f} kcal/mol")
        print(f"Range: [{min(benign_ddg):.2f}, {max(benign_ddg):.2f}] kcal/mol")

        effects = {}
        for r in benign:
            effects[r['effect']] = effects.get(r['effect'], 0) + 1
        print("Effect distribution:", effects)

    print("\n" + "-" * 70)
    print("STATISTICAL COMPARISON")
    print("-" * 70)
    if path_ddg and benign_ddg:
        mean_diff = statistics.mean(path_ddg) - statistics.mean(benign_ddg)
        print(f"Mean difference (Pathogenic - Benign): {mean_diff:.2f} kcal/mol")

        # Cohen's d
        pooled_std = ((statistics.stdev(path_ddg)**2 * (len(path_ddg)-1) +
                       statistics.stdev(benign_ddg)**2 * (len(benign_ddg)-1)) /
                      (len(path_ddg) + len(benign_ddg) - 2)) ** 0.5
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
        print(f"Effect size (Cohen's d): {cohens_d:.2f}")

        if cohens_d < 0.2:
            effect_interpretation = "negligible"
        elif cohens_d < 0.5:
            effect_interpretation = "small"
        elif cohens_d < 0.8:
            effect_interpretation = "medium"
        else:
            effect_interpretation = "large"
        print(f"Effect interpretation: {effect_interpretation}")

    # Threshold analysis
    print("\n" + "-" * 70)
    print("THRESHOLD ANALYSIS FOR PATHOGENICITY PREDICTION")
    print("-" * 70)

    thresholds = [0.5, 1.0, 1.5, 2.0, 2.5]
    for thresh in thresholds:
        path_above = sum(1 for d in path_ddg if d > thresh)
        benign_above = sum(1 for d in benign_ddg if d > thresh)
        sensitivity = path_above / len(path_ddg) if path_ddg else 0
        specificity = (len(benign_ddg) - benign_above) / len(benign_ddg) if benign_ddg else 0
        balanced_acc = (sensitivity + specificity) / 2

        print(f"\nThreshold: ddG > {thresh} kcal/mol")
        print(f"  Sensitivity: {sensitivity*100:.1f}% ({path_above}/{len(path_ddg)} pathogenic)")
        print(f"  Specificity: {specificity*100:.1f}% ({len(benign_ddg)-benign_above}/{len(benign_ddg)} benign)")
        print(f"  Balanced Accuracy: {balanced_acc*100:.1f}%")

    print("\n" + "=" * 70)
    print(f"Results saved to: {OUTPUT_FILE}")
    print("=" * 70)

if __name__ == "__main__":
    main()
