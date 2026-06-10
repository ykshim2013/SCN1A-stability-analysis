#!/usr/bin/env python3
"""
Create manuscript figures for SCN1A stability analysis paper (v7 — revised for SCI journal).

Generates three publication-quality figures + one supplementary figure:
- Figure 1: Stability analysis of 1,011 variants (violin + pLDDT ROC + correlation + ROC)
- Figure 2: Domain-stratified analysis (dual-axis box-strip + grouped bar enrichment)
- Figure 3: Nav1.1 linear domain architecture with variant stability landscape (standalone)
- Supplementary Figure S1: Predictor comparison ROC

Changes from v6:
  - Removed old Figure 3 (functional class analysis) — negative results conveyed by Table 1
  - Old Figure 4 → Figure 3, completely redesigned as professional linear architecture
  - Removed Figure 3 Panel B (redundant with Figure 2B) — Figure 3 is now standalone topology
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle
from sklearn.metrics import roc_curve, auc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ─── Publication-quality defaults ───
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['figure.dpi'] = 300

# ─── Data paths (repo-relative) ───
# Input calibration tables are NOT distributed in this repository prior to
# publication; place them under data/ to reproduce the figures locally.
import os
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE = os.path.join(_REPO, 'data', 'calibration_set.tsv')
EXT_PRED_FILE = os.path.join(_REPO, 'data', 'external_predictors.tsv')
OUTPUT_DIR = os.path.join(_REPO, 'figures') + os.sep
FIG_PREFIX = ''

# ─── Load datasets ───
df = pd.read_csv(DATA_FILE, sep='\t')
df['Label'] = (df['Classification'] == 'Pathogenic').astype(int)
print(f"Loaded {len(df)} variants (Pathogenic: {df['Label'].sum()}, Benign: {(df['Label']==0).sum()})")

df_ext = pd.read_csv(EXT_PRED_FILE, sep='\t')
print(f"Loaded external predictors: {len(df_ext)} variants")

# ─── Domain region assignment ───
def assign_region(row):
    seg_type = row['Segment_Type'] if 'Segment_Type' in row.index else row.get('Segment_type', '')
    if seg_type == 'VSD_S4':
        return 'S4 voltage sensor'
    if seg_type == 'VSD':
        return 'VSD helices (S1-S3)'
    if seg_type == 'Pore_helix':
        return 'Pore helices (S5-S6)'
    if seg_type == 'P_loop':
        return 'P-loops / select. filter'
    if seg_type == 'Loop/Terminus':
        return 'Loops / termini'
    return 'Other'

df['Region'] = df.apply(assign_region, axis=1)

REGION_ORDER = [
    'S4 voltage sensor',
    'VSD helices (S1-S3)',
    'Pore helices (S5-S6)',
    'P-loops / select. filter',
    'Loops / termini',
]
DOMAIN_COLORS = ['#FF9800', '#42A5F5', '#66BB6A', '#AB47BC', '#78909C']


def _add_sig_bracket(ax, x1, x2, y, p, h=0.5, fs=9):
    """Draw significance bracket between two x-positions."""
    if p < 0.001:
        sig = '***'
    elif p < 0.01:
        sig = '**'
    elif p < 0.05:
        sig = '*'
    else:
        sig = 'n.s.'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.0, color='black', clip_on=False)
    ax.text((x1+x2)/2, y+h+0.1, sig, ha='center', va='bottom', fontsize=fs,
            fontweight='bold', clip_on=False)


# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Protein stability analysis of 1,011 SCN1A variants
# ═══════════════════════════════════════════════════════════════
def create_figure1():
    """
    Figure 1: (A) pLDDT-stratified ROC, (B) Violin plots, (C) Correlation, (D) ROC curves
    """
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.4)
    axes = [fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[1, 1])]

    fx_path = df[df['Label'] == 1]['FoldX_ddG'].dropna()
    fx_ben = df[df['Label'] == 0]['FoldX_ddG'].dropna()
    ro_path = df[df['Label'] == 1]['Rosetta_ddG_REU'].dropna()
    ro_ben = df[df['Label'] == 0]['Rosetta_ddG_REU'].dropna()

    col_ben, col_path = '#2ecc71', '#e74c3c'

    # ── Panel B: Dual-axis violin plots ──
    ax = axes[1]
    fx_data = [fx_ben.values, fx_path.values]
    parts_fx = ax.violinplot(fx_data, positions=[0.8, 1.6], showmeans=False,
                             showmedians=False, showextrema=False, widths=0.55)
    for i, pc in enumerate(parts_fx['bodies']):
        pc.set_facecolor([col_ben, col_path][i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)

    bp_fx = ax.boxplot(fx_data, positions=[0.8, 1.6], widths=0.12,
                       patch_artist=True, showfliers=False)
    for patch in bp_fx['boxes']:
        patch.set_facecolor('white'); patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp_fx[element], color='black')

    ax.set_ylabel('FoldX ΔΔG (kcal/mol)', fontsize=11, color='#333333')
    ax.set_ylim(-10, 28)
    ax.spines['top'].set_visible(False)
    _, p_fx = stats.mannwhitneyu(fx_ben, fx_path)
    _add_sig_bracket(ax, 0.8, 1.6, 29, p_fx, h=0.5, fs=11)

    ax2 = ax.twinx()
    ro_data = [ro_ben.values, ro_path.values]
    parts_ro = ax2.violinplot(ro_data, positions=[2.8, 3.6], showmeans=False,
                              showmedians=False, showextrema=False, widths=0.55)
    for i, pc in enumerate(parts_ro['bodies']):
        pc.set_facecolor([col_ben, col_path][i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
        pc.set_linestyle('--')

    bp_ro = ax2.boxplot(ro_data, positions=[2.8, 3.6], widths=0.12,
                        patch_artist=True, showfliers=False)
    for patch in bp_ro['boxes']:
        patch.set_facecolor('white'); patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp_ro[element], color='black')

    ax2.set_ylabel('Rosetta Cartesian ΔΔG (REU)', fontsize=11, color='#333333')
    ax2.set_ylim(-20, 45)
    ax2.spines['top'].set_visible(False)
    _, p_ro = stats.mannwhitneyu(ro_ben, ro_path)
    _add_sig_bracket(ax2, 2.8, 3.6, 47, p_ro, h=0.8, fs=11)

    ax.set_xticks([0.8, 1.6, 2.8, 3.6])
    ax.set_xticklabels(['B/LB', 'P/LP', 'B/LB', 'P/LP'], fontsize=9)
    ax.set_xlim(0.2, 4.2)
    ax.text(-0.12, 1.05, 'B', fontsize=14, fontweight='bold', transform=ax.transAxes)
    ben_patch = mpatches.Patch(color=col_ben, alpha=0.6, label=f'B/LB (n={len(fx_ben)})')
    path_patch = mpatches.Patch(color=col_path, alpha=0.6, label=f'P/LP (n={len(fx_path)})')
    ax.legend(handles=[ben_patch, path_patch], loc='upper center', fontsize=8, framealpha=0.9)

    # ── Panel A: pLDDT-stratified ROC ──
    ax = axes[0]
    # R2: pLDDT already present in DATA_FILE, so use df directly
    df_merged = df.copy()
    df_merged = df_merged.dropna(subset=['pLDDT', 'FoldX_ddG', 'Rosetta_ddG_REU'])

    plddt_bins = [
        ('pLDDT < 70\n(low)', df_merged['pLDDT'] < 70),
        ('70 ≤ pLDDT < 90\n(confident)', (df_merged['pLDDT'] >= 70) & (df_merged['pLDDT'] < 90)),
        ('pLDDT ≥ 90\n(very high)', df_merged['pLDDT'] >= 90),
    ]

    bar_aucs_fx, bar_aucs_ro, bar_aucs_comb, bar_labels, bar_n = [], [], [], [], []
    for label, mask in plddt_bins:
        subset = df_merged[mask]
        n_path = subset['Label'].sum()
        n_ben = (subset['Label'] == 0).sum()
        if n_ben >= 2 and n_path >= 2:
            y_true = subset['Label'].values
            fx_s = subset['FoldX_ddG'].values
            ro_s = subset['Rosetta_ddG_REU'].values
            fx_n = (fx_s - fx_s.mean()) / (fx_s.std() + 1e-10)
            ro_n = (ro_s - ro_s.mean()) / (ro_s.std() + 1e-10)
            fpr, tpr, _ = roc_curve(y_true, fx_s); auc_fx = auc(fpr, tpr)
            fpr, tpr, _ = roc_curve(y_true, ro_s); auc_ro = auc(fpr, tpr)
            fpr, tpr, _ = roc_curve(y_true, (fx_n+ro_n)/2); auc_c = auc(fpr, tpr)
        else:
            auc_fx = auc_ro = auc_c = np.nan
        bar_aucs_fx.append(auc_fx); bar_aucs_ro.append(auc_ro)
        bar_aucs_comb.append(auc_c); bar_labels.append(label); bar_n.append(len(subset))

    x_pos = np.arange(len(bar_labels))
    width = 0.25
    ax.bar(x_pos - width, bar_aucs_fx, width, color='#3498db', alpha=0.8, edgecolor='white', label='FoldX')
    ax.bar(x_pos, bar_aucs_ro, width, color='#9b59b6', alpha=0.8, edgecolor='white', label='Rosetta')
    ax.bar(x_pos + width, bar_aucs_comb, width, color='#e67e22', alpha=0.8, edgecolor='white', label='Combined')
    for i, n in enumerate(bar_n):
        vals = [bar_aucs_fx[i], bar_aucs_ro[i], bar_aucs_comb[i]]
        max_v = max(v for v in vals if not np.isnan(v)) if any(not np.isnan(v) for v in vals) else 0.5
        ax.text(i, max_v + 0.03, f'n={n}', ha='center', va='bottom', fontsize=8, fontstyle='italic')
    ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='Random')
    ax.set_ylabel('ROC-AUC', fontsize=11)
    ax.set_xlabel('AlphaFold3 structural confidence', fontsize=10)
    ax.set_xticks(x_pos); ax.set_xticklabels(bar_labels, fontsize=7.5)
    ax.set_ylim(0, 1.1)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left', fontsize=8, framealpha=0.95)
    ax.text(-0.12, 1.05, 'A', fontsize=14, fontweight='bold', transform=ax.transAxes)

    # ── Panel C: Correlation scatter ──
    ax = axes[2]
    valid_scatter = df.dropna(subset=['FoldX_ddG', 'Rosetta_ddG_REU'])
    path_mask = valid_scatter['Label'] == 1
    ben_mask = valid_scatter['Label'] == 0
    ax.scatter(valid_scatter.loc[path_mask, 'FoldX_ddG'],
               valid_scatter.loc[path_mask, 'Rosetta_ddG_REU'],
               color='#e74c3c', alpha=0.3, s=15, edgecolors='none', label='P/LP')
    ax.scatter(valid_scatter.loc[ben_mask, 'FoldX_ddG'],
               valid_scatter.loc[ben_mask, 'Rosetta_ddG_REU'],
               color='#2ecc71', alpha=0.5, s=25, edgecolors='none', label='B/LB')
    ax.axhline(y=0, color='black', linewidth=0.3, alpha=0.3)
    ax.axvline(x=0, color='black', linewidth=0.3, alpha=0.3)
    ax.set_xlabel('FoldX ΔΔG (kcal/mol)', fontsize=11)
    ax.set_ylabel('Rosetta Cartesian ΔΔG (REU)', fontsize=11)
    ax.text(-0.12, 1.05, 'C', fontsize=14, fontweight='bold', transform=ax.transAxes)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)

    # ── Panel D: ROC curves ──
    ax = axes[3]
    valid_data = df.dropna(subset=['FoldX_ddG', 'Rosetta_ddG_REU', 'Label'])
    y_true = valid_data['Label'].values
    foldx_s = valid_data['FoldX_ddG'].values
    rosetta_s = valid_data['Rosetta_ddG_REU'].values
    fx_n = (foldx_s - foldx_s.mean()) / foldx_s.std()
    ro_n = (rosetta_s - rosetta_s.mean()) / rosetta_s.std()
    comb_s = (fx_n + ro_n) / 2

    for scores, color, label in [(foldx_s, '#3498db', 'FoldX'),
                                  (rosetta_s, '#9b59b6', 'Rosetta'),
                                  (comb_s, '#e67e22', 'Combined')]:
        fpr, tpr, _ = roc_curve(y_true, scores)
        a = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=color, linewidth=2.5, label=f'{label} (AUC = {a:.3f})')
        j = tpr - fpr
        opt = np.argmax(j)
        ax.scatter(fpr[opt], tpr[opt], color=color, s=80, zorder=5, edgecolors='black', linewidth=1.5)

    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
    ax.set_xlabel('False Positive Rate (1 − Specificity)', fontsize=11)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=11)
    ax.text(-0.12, 1.05, 'D', fontsize=14, fontweight='bold', transform=ax.transAxes)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.95)
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(f'{OUTPUT_DIR}2606_Figure1_stability_analysis.{ext}', dpi=300,
                    bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(f'{OUTPUT_DIR}2606_Figure1_stability_analysis.tiff', dpi=600,
                bbox_inches='tight', facecolor='white', edgecolor='none',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print("✓ Figure 1: Stability analysis")


# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Domain-stratified prediction performance
# ═══════════════════════════════════════════════════════════════
def create_figure2():
    """
    Figure 2: Dual-axis box-strip of FoldX (left) and Rosetta (right) ΔΔG
    distributions for P/LP variants, stratified by Nav1.1 functional domain.
    The previous Panel B (stacked composition bars) was removed in revision;
    the same within-domain enrichment information is reported in Table 1 and
    Supplementary Table S3.
    """
    fig, ax = plt.subplots(1, 1, figsize=(9, 6.5))
    pathogenic_df = df[df['Label'] == 1].copy()

    # ── Dual-axis box-strip (formerly Panel A) ──
    n_regions = len(REGION_ORDER)
    pos_fx = np.arange(n_regions) * 2
    pos_ro = np.arange(n_regions) * 2 + 0.7

    for i, (region, color) in enumerate(zip(REGION_ORDER, DOMAIN_COLORS)):
        vals = pathogenic_df[pathogenic_df['Region'] == region]['FoldX_ddG'].dropna().values
        if len(vals) == 0: continue
        bp = ax.boxplot([vals], positions=[pos_fx[i]], widths=0.5, patch_artist=True,
                        showfliers=False, zorder=2)
        bp['boxes'][0].set_facecolor(color); bp['boxes'][0].set_alpha(0.35)
        bp['boxes'][0].set_edgecolor(color)
        bp['medians'][0].set_color(color); bp['medians'][0].set_linewidth(2)
        for w in bp['whiskers']: w.set_color(color)
        for c in bp['caps']: c.set_color(color)
        rng = np.random.default_rng(42 + i)
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(pos_fx[i] + jitter, vals, color=color, edgecolor='white',
                   s=10, zorder=3, linewidth=0.3, alpha=0.4)

    ax.set_ylabel('FoldX ΔΔG (kcal/mol)', fontsize=11)
    ax.set_ylim(-8, 22)
    ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.3)
    ax.spines['top'].set_visible(False)

    ax2 = ax.twinx()
    for i, (region, color) in enumerate(zip(REGION_ORDER, DOMAIN_COLORS)):
        vals = pathogenic_df[pathogenic_df['Region'] == region]['Rosetta_ddG_REU'].dropna().values
        if len(vals) == 0: continue
        bp = ax2.boxplot([vals], positions=[pos_ro[i]], widths=0.5, patch_artist=True,
                         showfliers=False, zorder=2)
        bp['boxes'][0].set_facecolor(color); bp['boxes'][0].set_alpha(0.15)
        bp['boxes'][0].set_edgecolor(color); bp['boxes'][0].set_linestyle('--')
        bp['medians'][0].set_color(color); bp['medians'][0].set_linewidth(2)
        bp['medians'][0].set_linestyle('--')
        for w in bp['whiskers']: w.set_color(color); w.set_linestyle('--')
        for c in bp['caps']: c.set_color(color); c.set_linestyle('--')
        rng = np.random.default_rng(142 + i)
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax2.scatter(pos_ro[i] + jitter, vals, color=color, marker='s',
                    edgecolor='white', s=10, zorder=3, linewidth=0.3, alpha=0.3)

    ax2.set_ylabel('Rosetta Cartesian ΔΔG (REU)', fontsize=11)
    ax2.set_ylim(-15, 45); ax2.spines['top'].set_visible(False)
    tick_positions = (pos_fx + pos_ro) / 2
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([r.replace(' / ', '/\n') for r in REGION_ORDER], fontsize=8)
    ax.set_xlim(-0.7, pos_ro[-1] + 0.7)

    s4_fx = pathogenic_df[pathogenic_df['Region'] == 'S4 voltage sensor']['FoldX_ddG'].dropna().values
    pl_fx = pathogenic_df[pathogenic_df['Region'] == 'P-loops / select. filter']['FoldX_ddG'].dropna().values
    _, pw_p = stats.mannwhitneyu(s4_fx, pl_fx, alternative='two-sided')
    if pw_p < 0.05:
        _add_sig_bracket(ax, pos_fx[0], pos_fx[3], 20, pw_p, h=1.0, fs=10)

    leg_fx = Line2D([0], [0], color='gray', linewidth=2, linestyle='-', label='FoldX')
    leg_ro = Line2D([0], [0], color='gray', linewidth=2, linestyle='--', label='Rosetta')
    ax.legend(handles=[leg_fx, leg_ro], loc='upper right', fontsize=8, framealpha=0.95)

    # Former Panel B (100%-stacked composition bars) was removed in
    # revision. The same within-domain enrichment information is
    # reported in Table 1 and Supplementary Table S3.

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(f'{OUTPUT_DIR}2606_Figure2_domain_analysis.{ext}', dpi=300,
                    bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(f'{OUTPUT_DIR}2606_Figure2_domain_analysis.tiff', dpi=600,
                bbox_inches='tight', facecolor='white', edgecolor='none',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print("✓ Figure 2: Domain-stratified analysis (revision: bar percentages displayed on figure)")


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Nav1.1 Linear Domain Architecture (PROFESSIONAL REDESIGN)
# ═══════════════════════════════════════════════════════════════
def create_figure3():
    """
    Figure 3: Nav1.1 variant stability landscape (standalone).

    Top: Variant lollipop plot (FoldX ΔΔG vs position) with domain brackets.
    Middle: Linear domain architecture bar colored by segment type.
    Bottom: Segment-level heatmap of mean ΔΔG.
    """
    pathogenic_df = df[df['Label'] == 1].copy()
    benign_df = df[df['Label'] == 0].copy()

    # ─── Build segment map with positions ───
    segments_info = []
    for seg_name in df['Segment'].unique():
        sub = df[df['Segment'] == seg_name]
        pos_min = sub['Position'].min()
        pos_max = sub['Position'].max()
        seg_type = sub['Segment_Type'].iloc[0]
        domain = sub['Domain'].iloc[0]
        n_path = (sub['Classification'] == 'Pathogenic').sum()
        n_ben = (sub['Classification'] != 'Pathogenic').sum()
        mean_ddg = sub[sub['Classification'] == 'Pathogenic']['FoldX_ddG'].mean() if n_path > 0 else 0
        segments_info.append({
            'segment': seg_name, 'pos_min': pos_min, 'pos_max': pos_max,
            'seg_type': seg_type, 'domain': domain,
            'n_path': n_path, 'n_ben': n_ben, 'mean_ddg': mean_ddg
        })
    seg_df = pd.DataFrame(segments_info).sort_values('pos_min').reset_index(drop=True)

    # ─── Color scheme by segment type ───
    seg_type_colors = {
        'VSD_S4': '#FF9800',      # Orange - voltage sensor
        'VSD': '#42A5F5',          # Blue - VSD helices
        'Pore_helix': '#66BB6A',   # Green - pore helices
        'P_loop': '#AB47BC',       # Purple - P-loops
        'Loop/Terminus': '#BDBDBD' # Gray - loops/termini
    }
    seg_type_labels = {
        'VSD_S4': 'S4 voltage sensor',
        'VSD': 'VSD helices (S1–S3)',
        'Pore_helix': 'Pore helices (S5–S6)',
        'P_loop': 'P-loop / selectivity filter',
        'Loop/Terminus': 'Loops / termini'
    }

    # R2: single-hue red sequential colormap for the bottom Mean ΔΔG heatmap
    # (per user request — all-red gradient, white → deep red).
    cmap = LinearSegmentedColormap.from_list('stability_reds',
        ['#FFFFFF', '#FFEBEE', '#FFCDD2', '#EF9A9A', '#E57373', '#EF5350',
         '#E53935', '#C62828', '#B71C1C'], N=256)
    norm = Normalize(vmin=0, vmax=8)

    # Domain boundaries for bracket annotations
    domain_ranges = {
        'DI': (132, 416), 'DII': (724, 976),
        'DIII': (1204, 1454), 'DIV': (1526, 1771)
    }

    protein_length = 2009

    # ═══ CREATE FIGURE (3 rows only — no Panel B) ═══
    # Revision: increased figure dimensions for label legibility per Reviewer 1
    # R2: domain architecture bar removed per user request — figure is now
    # lollipop (top) + mean-ΔΔG heatmap (bottom).
    fig = plt.figure(figsize=(18, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.0, 0.7],
                          hspace=0.10, left=0.07, right=0.84, top=0.94, bottom=0.10)

    # ── Variant lollipop plot (FoldX ΔΔG vs position) ──
    ax_lollipop = fig.add_subplot(gs[0])

    # Plot all pathogenic variants as lollipops
    for _, row in pathogenic_df.iterrows():
        pos = row['Position']
        ddg = row['FoldX_ddG']
        seg_type = row['Segment_Type']
        color = seg_type_colors.get(seg_type, '#999999')
        ax_lollipop.plot([pos, pos], [0, ddg], color=color, linewidth=0.4, alpha=0.5, zorder=1)
        ax_lollipop.scatter(pos, ddg, color=color, s=8, alpha=0.6, edgecolors='none', zorder=2)

    # Benign variants as green diamonds (larger for visibility)
    for _, row in benign_df.iterrows():
        pos = row['Position']
        ddg = row['FoldX_ddG']
        ax_lollipop.scatter(pos, ddg, color='#2ecc71', s=35, marker='D',
                           alpha=0.85, edgecolors='white', linewidth=0.6, zorder=3)

    # ── Annotate key variants discussed in manuscript text ──
    y_max = 10  # Must match the y-axis clip value used below
    key_variants = {
        'R1648H': {'ha': 'left', 'va': 'bottom', 'dx': 20, 'dy': 0.4},    # GOF, high ΔΔG
        'R101W':  {'ha': 'left', 'va': 'top', 'dx': 20, 'dy': -0.5},    # cLOF, clipped — annotate below bracket
        'T782I':  {'ha': 'right', 'va': 'top', 'dx': -20, 'dy': -0.5},    # GOF, stabilizing extreme
        'R1639C': {'ha': 'left', 'va': 'bottom', 'dx': 20, 'dy': 0.5},    # S4 gating charge, near-neutral
        'R1645Q': {'ha': 'right', 'va': 'bottom', 'dx': -20, 'dy': 0.5},  # S4 gating charge, near-neutral
        'K1313T': {'ha': 'left', 'va': 'bottom', 'dx': 20, 'dy': 0.5},    # S4 gating charge, near-neutral
    }
    all_df = pd.concat([pathogenic_df, benign_df])
    for var_name, offset in key_variants.items():
        match = all_df[all_df['Variant'] == var_name]
        if match.empty:
            continue
        r = match.iloc[0]
        ddg = r['FoldX_ddG']
        # For clipped outliers, annotate at the clipped position with ddG value
        if ddg > y_max:
            label = f"{var_name} (ΔΔG {ddg:.1f})"
            xy_y = y_max - 0.15
        else:
            label = var_name
            xy_y = ddg
        ax_lollipop.annotate(
            label, xy=(r['Position'], xy_y),
            xytext=(r['Position'] + offset['dx'], xy_y + offset['dy']),
            fontsize=13, fontstyle='italic', color='#222222', fontweight='bold',
            ha=offset['ha'], va=offset['va'],
            arrowprops=dict(arrowstyle='-', color='#444444', linewidth=0.8, shrinkA=0, shrinkB=2),
            zorder=5
        )

    # Reference lines (R2 revision: red line now marks the Youden-optimal
    # discrimination threshold (ΔΔG = 1.60 kcal/mol), not the 2.0 strongly-
    # destabilising cutoff; the ±0.5 near-neutral grey lines were removed.)
    ax_lollipop.axhline(y=0, color='black', linewidth=1.0, alpha=0.7)
    ax_lollipop.axhline(y=1.60, color='#D32F2F', linewidth=1.5, linestyle='-', alpha=0.55)

    # Clip extreme outliers visually — place triangle markers at top edge
    y_max = 10
    clipped = pathogenic_df[pathogenic_df['FoldX_ddG'] > y_max]
    for _, row in clipped.iterrows():
        ax_lollipop.scatter(row['Position'], y_max - 0.15, marker='^', color='#B71C1C',
                           s=18, zorder=4, edgecolors='none', alpha=0.8)

    ax_lollipop.set_ylabel('FoldX ΔΔG\n(kcal/mol)', fontsize=15, fontweight='bold')
    ax_lollipop.tick_params(axis='both', labelsize=12)
    ax_lollipop.tick_params(axis='y', labelsize=10)
    ax_lollipop.set_xlim(0, protein_length + 10)
    ax_lollipop.set_ylim(-4, y_max)
    ax_lollipop.spines['top'].set_visible(False)
    ax_lollipop.spines['right'].set_visible(False)
    ax_lollipop.spines['bottom'].set_visible(False)
    ax_lollipop.set_xticks([])

    # Domain brackets at top (above clipped y-axis, using clip_on=False)
    for dom_label, (start, end) in domain_ranges.items():
        mid = (start + end) / 2
        y_bracket = y_max + 0.3
        ax_lollipop.plot([start, start, end, end],
                        [y_bracket - 0.2, y_bracket, y_bracket, y_bracket - 0.2],
                        color='#333333', linewidth=1.5, clip_on=False)
        ax_lollipop.text(mid, y_bracket + 0.2, dom_label, ha='center', va='bottom',
                        fontsize=16, fontweight='bold', color='#222222', clip_on=False)

    # N and C terminus labels
    ax_lollipop.text(-15, -1, 'N', fontsize=15, fontweight='bold', color='#1976D2',
                     ha='center', va='center')
    ax_lollipop.text(protein_length + 15, -1, 'C', fontsize=15, fontweight='bold',
                     color='#D32F2F', ha='center', va='center')

    # Linker region shading (strengthened for visibility)
    linker_regions = [(420, 707), (977, 1196), (1454, 1526)]
    for start, end in linker_regions:
        ax_lollipop.axvspan(start, end, alpha=0.12, color='#B0BEC5', zorder=0)

    # ── Lollipop legend (variant markers + reference lines) — right of top panel ──
    # R2: legend box and markers enlarged per user request.
    lollipop_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=seg_type_colors['VSD_S4'],
               markersize=11, label=seg_type_labels['VSD_S4']),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=seg_type_colors['VSD'],
               markersize=11, label=seg_type_labels['VSD']),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=seg_type_colors['Pore_helix'],
               markersize=11, label=seg_type_labels['Pore_helix']),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=seg_type_colors['P_loop'],
               markersize=11, label=seg_type_labels['P_loop']),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=seg_type_colors['Loop/Terminus'],
               markersize=11, label=seg_type_labels['Loop/Terminus']),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#2ecc71',
               markersize=12, label='B/LB variants'),
        Line2D([0], [0], color='#D32F2F', linewidth=2.4, linestyle='-', alpha=0.7,
               label='Youden-optimal threshold (ΔΔG = 1.60)'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='#B71C1C',
               markersize=11, label='Clipped higher ΔΔG (> 10 kcal/mol)'),
    ]
    ax_lollipop.legend(handles=lollipop_handles, loc='upper left', fontsize=14,
                       framealpha=0.95, ncol=1, handletextpad=0.7,
                       labelspacing=0.85, borderpad=0.8,
                       bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
                       edgecolor='#888888')

    # ── Stability heatmap bar (mean ΔΔG per segment) ──
    # R2: domain architecture bar removed per user request; heatmap promoted
    # from gs[2] → gs[1].
    ax_heat = fig.add_subplot(gs[1], sharex=ax_lollipop)

    for _, row in seg_df.iterrows():
        width = row['pos_max'] - row['pos_min'] + 1
        color = cmap(norm(row['mean_ddg'])) if row['n_path'] > 0 else '#F5F5F5'
        rect = Rectangle((row['pos_min'], 0), width, 1,
                         facecolor=color, edgecolor='white', linewidth=0.3, zorder=2)
        ax_heat.add_patch(rect)

    ax_heat.set_ylim(0, 1)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel('Mean\nΔΔG', fontsize=13, rotation=0, labelpad=35, va='center', fontweight='bold')
    ax_heat.spines['top'].set_visible(False)
    ax_heat.spines['right'].set_visible(False)
    ax_heat.spines['left'].set_visible(False)
    ax_heat.set_xlabel('Amino acid position', fontsize=15, fontweight='bold')
    ax_heat.tick_params(axis='x', labelsize=10)

    # ── Colorbar for heatmap — right of bottom panel ──
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    # R2: colorbar height matches the heatmap panel height (user request).
    # Pull the heatmap axis bbox after layout to anchor the colorbar to its
    # actual rendered geometry.
    fig.canvas.draw()
    pos_heat = ax_heat.get_position()
    cbar_ax = fig.add_axes([0.86, pos_heat.y0, 0.016, pos_heat.height])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Mean ΔΔG\n(kcal/mol)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)

    for ext in ['png', 'pdf']:
        plt.savefig(f'{OUTPUT_DIR}2606_Supplementary_Figure_S1_topology.{ext}', dpi=300,
                    bbox_inches='tight', facecolor='white', edgecolor='none')
    # 600 DPI TIFF for journal submission
    plt.savefig(f'{OUTPUT_DIR}2606_Supplementary_Figure_S1_topology.tiff', dpi=600,
                bbox_inches='tight', facecolor='white', edgecolor='none',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print("✓ Figure 3: Nav1.1 variant stability landscape (revision: larger fonts, solid reference lines)")


# ═══════════════════════════════════════════════════════════════
# SUPPLEMENTARY FIGURE S1: Predictor comparison ROC
# ═══════════════════════════════════════════════════════════════
def create_supp_figure_s1():
    """Supplementary Figure S1: ROC comparison of all predictors."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    df_m = df.merge(df_ext[['Variant', 'am_pathogenicity', 'EVE_score']].drop_duplicates(),
                    on='Variant', how='left')
    df_m = df_m.dropna(subset=['FoldX_ddG', 'Rosetta_ddG_REU', 'Label'])
    y_true = df_m['Label'].values

    ax = axes[0]
    predictors = []

    fx = df_m['FoldX_ddG'].values
    fpr, tpr, _ = roc_curve(y_true, fx); a = auc(fpr, tpr)
    predictors.append(('FoldX', fpr, tpr, a, '#3498db', '-'))

    ro = df_m['Rosetta_ddG_REU'].values
    fpr, tpr, _ = roc_curve(y_true, ro); a = auc(fpr, tpr)
    predictors.append(('Rosetta', fpr, tpr, a, '#9b59b6', '-'))

    fx_n = (fx - fx.mean()) / fx.std()
    ro_n = (ro - ro.mean()) / ro.std()
    fpr, tpr, _ = roc_curve(y_true, (fx_n+ro_n)/2); a = auc(fpr, tpr)
    predictors.append(('Combined stability', fpr, tpr, a, '#e67e22', '-'))

    am_valid = df_m.dropna(subset=['am_pathogenicity'])
    if len(am_valid) > 10:
        fpr, tpr, _ = roc_curve(am_valid['Label'].values, am_valid['am_pathogenicity'].values)
        a = auc(fpr, tpr)
        predictors.append(('AlphaMissense', fpr, tpr, a, '#e74c3c', '--'))

    eve_valid = df_m.dropna(subset=['EVE_score'])
    if len(eve_valid) > 10:
        fpr, tpr, _ = roc_curve(eve_valid['Label'].values, -eve_valid['EVE_score'].values)
        a = auc(fpr, tpr)
        predictors.append(('EVE', fpr, tpr, a, '#2ecc71', '--'))

    predictors.sort(key=lambda x: x[3], reverse=True)
    for name, fpr, tpr, a, color, ls in predictors:
        ax.plot(fpr, tpr, color=color, linewidth=2.5, linestyle=ls, label=f'{name} (AUC = {a:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
    ax.set_xlabel('False Positive Rate (1 − Specificity)', fontsize=11)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=11)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.95)
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_title('Pathogenicity discrimination: all predictors', fontsize=11, fontweight='bold')
    ax.text(-0.12, 1.05, 'A', fontsize=14, fontweight='bold', transform=ax.transAxes)

    ax = axes[1]
    pred_names = [p[0] for p in predictors]
    pred_aucs = [p[3] for p in predictors]
    pred_colors = [p[4] for p in predictors]
    ax.barh(range(len(pred_names)), pred_aucs, color=pred_colors, alpha=0.8,
            edgecolor='white', height=0.6)
    ax.set_yticks(range(len(pred_names)))
    ax.set_yticklabels(pred_names, fontsize=10)
    ax.set_xlabel('ROC-AUC', fontsize=11); ax.set_xlim(0.5, 1.0)
    ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    for i, a in enumerate(pred_aucs):
        ax.text(a + 0.005, i, f'{a:.3f}', va='center', fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_title('Predictor AUC comparison', fontsize=11, fontweight='bold')
    ax.text(-0.15, 1.05, 'B', fontsize=14, fontweight='bold', transform=ax.transAxes)
    ax.invert_yaxis()

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(f'{OUTPUT_DIR}2606_Supplementary_Figure_S1_predictor_comparison.{ext}', dpi=300,
                    bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Supplementary Figure S1: Predictor comparison")


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("=" * 70)
    print("Creating manuscript figures for SCN1A stability analysis (v7)")
    print("=" * 70)

    create_figure1()
    create_figure2()
    create_figure3()
    create_supp_figure_s1()

    print("\n" + "=" * 70)
    print("All figures saved to:", OUTPUT_DIR)
    print("=" * 70)
