"""
Generate the new supplementary figures and tables for the revision:

  Supplementary Figure S2 — gnomAD AF-stratified FoldX ΔΔG distribution
                            overlaid with ClinVar B/LB and P/LP.
  Supplementary Figure S3 — ClinVar VUS FoldX ΔΔG distribution
                            stratified by AlphaMissense classification.
  Supplementary Table S2 (expanded) — full per-variant master table
                                      (P/LP + B/LB + VUS + gnomAD).
  Supplementary Table S6 — gnomAD AF-stratified summary statistics.
  Supplementary Table S7 — cryo-EM Nav structural-state comparison.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
FX = ROOT / "foldx_analysis"
RESULTS = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SUPPL = ROOT
DATA_PATH = FX / "master_expanded_dataset.tsv"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 1.0,
    "savefig.dpi": 300,
})


def save_three(fig, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


# ----------------------------------------------------------------------
# Supplementary Figure S2 — gnomAD distribution
# ----------------------------------------------------------------------
def make_s2(df: pd.DataFrame) -> None:
    plot_strata = [
        ("PLP",                "ClinVar P/LP",          "#D32F2F"),
        ("BLB",                "ClinVar B/LB",          "#2E7D32"),
        ("gnomAD_common",      "gnomAD common\n" + r"(AF $\geq$ 1 $\times$ 10$^{-4}$)",                                 "#1565C0"),
        ("gnomAD_rare",        "gnomAD rare\n" + r"(1 $\times$ 10$^{-5}$ $\leq$ AF $<$ 1 $\times$ 10$^{-4}$)",         "#42A5F5"),
        ("gnomAD_ultra_rare",  "gnomAD ultra-rare\n" + r"(AF $<$ 1 $\times$ 10$^{-5}$)",                               "#90CAF9"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    panel_letters = ["A", "B"]
    for ax_idx, (ax, metric, ylabel) in enumerate([
        (axes[0], "FoldX_ddG", "FoldX ΔΔG (kcal/mol)"),
        (axes[1], "am_pathogenicity", "AlphaMissense pathogenicity"),
    ]):
        data, labels, colors = [], [], []
        for stratum, label, color in plot_strata:
            sub = df[(df["Stratum"] == stratum) & df[metric].notna()][metric].values
            if len(sub) >= 5:
                data.append(sub)
                labels.append(f"{label}\n(n={len(sub)})")
                colors.append(color)
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for pc, c in zip(parts["bodies"], colors):
            pc.set_facecolor(c)
            pc.set_alpha(0.55)
            pc.set_edgecolor(c)
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            if key in parts:
                parts[key].set_edgecolor("#444444")
                parts[key].set_linewidth(1.0)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_ylabel(ylabel)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if metric == "FoldX_ddG":
            ax.axhline(0, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
            ax.axhline(2.0, color="#D32F2F", linewidth=1.0, alpha=0.5,
                       label="strongly destabilising (2.0 kcal/mol)")
            ax.legend(fontsize=8, loc="upper right")
        ax.text(-0.12, 1.05, panel_letters[ax_idx], fontsize=14,
                fontweight="bold", transform=ax.transAxes)
    fig.tight_layout()
    save_three(fig, FIG / "2606_Figure5_gnomad_AF_tier")
    plt.close(fig)
    print(f"[S1] wrote {FIG / '2606_Figure5_gnomad_AF_tier.{png,pdf,tiff}'}")


# ----------------------------------------------------------------------
# Supplementary Figure S3 — VUS distribution + AM overlay
# ----------------------------------------------------------------------
def make_s3(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    categories = [
        ("BLB", "B/LB", "#2E7D32"),
        ("VUS", "VUS",   "#FFA726"),
        ("PLP", "P/LP",  "#D32F2F"),
    ]

    # (A) Violin of FoldX ΔΔG for B/LB, VUS, P/LP
    data, labels, colors = [], [], []
    counts = {}
    for stratum, label, color in categories:
        sub = df[(df["Stratum"] == stratum) & df["FoldX_ddG"].notna()]["FoldX_ddG"].values
        data.append(sub)
        labels.append(f"{label}\n(n={len(sub)})")
        colors.append(color)
        counts[stratum] = len(sub)
    parts = axes[0].violinplot(data, showmedians=True, widths=0.75)
    for pc, c in zip(parts["bodies"], colors):
        pc.set_facecolor(c)
        pc.set_alpha(0.55)
        pc.set_edgecolor(c)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        if key in parts:
            parts[key].set_color("#444444")
            parts[key].set_linewidth(1.0)
    axes[0].set_xticks(range(1, len(labels) + 1))
    axes[0].set_xticklabels(labels, fontsize=10)
    axes[0].set_ylabel("FoldX ΔΔG (kcal/mol)")
    axes[0].axhline(0, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    axes[0].axhline(2.0, color="#D32F2F", linewidth=1.0, alpha=0.5,
                    label="strongly destabilising (2.0 kcal/mol)")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)
    axes[0].text(-0.12, 1.05, "A", fontsize=14, fontweight="bold",
                 transform=axes[0].transAxes)

    # (B) Within-category composition: stable / mild / strongly destabilising
    bins = [
        ("stabilising / near-neutral (ΔΔG ≤ 0.5)", "#2E7D32", lambda v: v <= 0.5),
        ("mildly destabilising (0.5 < ΔΔG ≤ 2.0)", "#FFA726", lambda v: (v > 0.5) & (v <= 2.0)),
        ("strongly destabilising (ΔΔG > 2.0)",     "#D32F2F", lambda v: v > 2.0),
    ]
    x_positions = np.arange(len(categories))
    bottom = np.zeros(len(categories))
    for bin_label, bin_color, bin_fn in bins:
        fractions = []
        for stratum, _, _ in categories:
            vals = df[(df["Stratum"] == stratum) & df["FoldX_ddG"].notna()]["FoldX_ddG"].values
            n = len(vals)
            frac = 100.0 * np.sum(bin_fn(vals)) / n if n else 0.0
            fractions.append(frac)
        bars = axes[1].bar(x_positions, fractions, bottom=bottom, width=0.6,
                           color=bin_color, edgecolor="white", linewidth=1.0,
                           label=bin_label)
        for bar, frac in zip(bars, fractions):
            if frac >= 4:
                axes[1].text(bar.get_x() + bar.get_width() / 2,
                             bar.get_y() + bar.get_height() / 2,
                             f"{frac:.1f}%", ha="center", va="center",
                             fontsize=9, color="white", fontweight="bold")
        bottom += np.array(fractions)
    axes[1].set_xticks(x_positions)
    axes[1].set_xticklabels([f"{label}\n(n={counts[stratum]})"
                             for stratum, label, _ in categories], fontsize=10)
    axes[1].set_ylabel("Within-category proportion (%)")
    axes[1].set_ylim(0, 100)
    axes[1].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)
    axes[1].text(-0.12, 1.05, "B", fontsize=14, fontweight="bold",
                 transform=axes[1].transAxes)

    fig.tight_layout()
    save_three(fig, FIG / "2606_Figure4_VUS_distribution")
    plt.close(fig)
    print(f"[S2] wrote {FIG / '2606_Figure4_VUS_distribution.{png,pdf,tiff}'}")


# ----------------------------------------------------------------------
# Supplementary Table S2 expanded — XLSX
# ----------------------------------------------------------------------
def make_table_s2(df: pd.DataFrame) -> None:
    out_xlsx = SUPPL / "2605_Supplementary_Table_S2_expanded.xlsx"
    cols = [
        "ProteinChange1", "Position", "Stratum", "ClinVar_class",
        "ReviewStatus", "ReviewCategory",
        "gnomAD_AC", "gnomAD_AN", "gnomAD_AF", "gnomAD_popmax_AF",
        "gnomAD_AF_tier", "Domain", "Segment", "Region",
        "FoldX_ddG", "FoldX_std", "FoldX_n_runs", "FoldX_category",
        "Rosetta_ddG", "am_pathogenicity", "am_class",
        "EVE_score", "EVE_class",
    ]
    cols = [c for c in cols if c in df.columns]
    sub = df[cols].copy()
    sub.to_excel(out_xlsx, index=False, sheet_name="all_variants")
    print(f"[S2 table] wrote {out_xlsx} ({len(sub)} variants × {len(cols)} cols)")


# ----------------------------------------------------------------------
# Supplementary Table S6 — gnomAD AF-stratified summary
# ----------------------------------------------------------------------
def make_table_s6(df: pd.DataFrame) -> None:
    rows = []
    for tier in ("ultra_rare", "rare", "common"):
        sub = df[df["Stratum"] == f"gnomAD_{tier}"]
        for col in ("FoldX_ddG", "am_pathogenicity"):
            v = sub[col].dropna()
            rows.append({
                "AF_tier": tier,
                "metric": col,
                "n": len(v),
                "mean": float(v.mean()) if len(v) else float("nan"),
                "median": float(v.median()) if len(v) else float("nan"),
                "std": float(v.std()) if len(v) > 1 else float("nan"),
                "fraction_FoldX_>2.0": float((sub["FoldX_ddG"] > 2.0).mean())
                                       if col == "FoldX_ddG" and len(v) else float("nan"),
            })
    summary = pd.DataFrame(rows)
    out_xlsx = SUPPL / "2605_Supplementary_Table_S6_gnomad.xlsx"
    with pd.ExcelWriter(out_xlsx) as writer:
        summary.to_excel(writer, sheet_name="AF_tier_summary", index=False)
        gn = df[df["Stratum"].str.startswith("gnomAD_", na=False)].copy()
        keep = ["ProteinChange1", "Position", "Stratum", "Domain", "Region",
                "gnomAD_AC", "gnomAD_AN", "gnomAD_AF", "gnomAD_popmax_AF",
                "FoldX_ddG", "FoldX_std", "am_pathogenicity", "am_class"]
        gn[[c for c in keep if c in gn.columns]].to_excel(
            writer, sheet_name="gnomAD_variants", index=False)
    print(f"[S6 table] wrote {out_xlsx}")


# ----------------------------------------------------------------------
# Supplementary Table S7 — cryo-EM Nav state comparison
# ----------------------------------------------------------------------
def make_table_s7() -> None:
    diag = pd.read_csv(RESULTS / "structural_state_diagnostic.tsv", sep="\t")
    out_xlsx = SUPPL / "2605_Supplementary_Table_S7_cryoem_comparison.xlsx"
    diag.to_excel(out_xlsx, index=False, sheet_name="state_diagnostic")
    print(f"[S7 table] wrote {out_xlsx}")


def main() -> None:
    if not DATA_PATH.exists():
        print(f"ERROR: {DATA_PATH} not yet built. Run integrate_expanded_dataset.py first.",
              file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(DATA_PATH, sep="\t")
    print(f"[load] {len(df)} variants from {DATA_PATH}")
    make_s2(df)
    make_s3(df)
    make_table_s2(df)
    make_table_s6(df)
    make_table_s7()


if __name__ == "__main__":
    main()
