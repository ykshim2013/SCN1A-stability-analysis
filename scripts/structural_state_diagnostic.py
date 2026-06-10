"""
Structural-state diagnostic of the SCN1A AlphaFold3 model used for
FoldX ΔΔG calculations.

Computes four metrics on `foldx_analysis/scn1a_wt_Repair.pdb` and on
deposited cryo-EM Nav structures, all interpretable in the same
geometric frame:

  (A) S4-helix gating-charge z-position per VSD (DI–DIV) relative
      to the DEKA selectivity-filter mid-plane. "Up" (activated) =
      large positive z; "down" (resting) ≈ 0 or negative.
  (B) IFM inactivation gate distance from receptor pocket: distance
      from F1489 Cα to the centroid of receptor residues
      (Y1494 + N1466 + N1494) — short = docked = inactivated.
  (C) Activation-gate pore radius: minimum perpendicular distance
      from the central pore axis to any heavy atom across S6
      bundle-crossing residues, scanned over a 30 Å z-range.
  (D) DEKA selectivity-filter Cα cross-section (D/E/K/A in DI/DII/
      DIII/DIV). Diagonal Cα–Cα distances diagnose pore patency.

Structures compared:
  - SCN1A AF3 (our model)
  - 7K48  NavAb/Nav1.7-VS2A trapped resting   [reference: resting]
  - 6J8E  Nav1.2-β2-KIIIA (pore-blocked)       [activated/inactivated]
  - 7XVE  Nav1.7 mutant class-I               [activated]
  - 6UZ3  Nav1.5 cardiac                       [inactivated]
  - 7DTC  Nav1.5 E1784K                        [inactivated mutant]
  - 7XM9  Nav1.7-β1-β2 + XEN907                [inactivated]
  - 6AGF  Nav1.4-β1                            [inactivated]

Outputs:
  - results/structural_state_diagnostic.tsv        wide table
  - results/pore_radius_profile.tsv                long form, per z
  - figures/Supplementary_Figure_S4_state_diagnostic.{png,pdf,tiff}
"""
from __future__ import annotations

import math
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB.Polypeptide import three_to_index, index_to_one

ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT / "structure" / "cryoem_nav_references"
AF3_PDB = ROOT / "foldx_analysis" / "scn1a_wt_Repair.pdb"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# SCN1A (Nav1.1, UniProt P35498) functional residues used as our reference frame.
# Other Nav isoforms have shifted numbering — see Pan/Yan-cohort residue maps
# below for matched residues.
SCN1A_DEKA = {  # Selectivity filter residues (DI, DII, DIII, DIV)
    "DI":   ("D", 382),
    "DII":  ("E", 951),
    "DIII": ("K", 1432),
    "DIV":  ("A", 1724),
}

# S4 helix arginine positions per VSD (verified against AF3 PDB).
# We compute the z-projection for ALL gating arginines and report
# the most extracellular (highest projected z) as R1.
SCN1A_S4_ARGS = {
    "DI":   [219, 222, 225, 228, 231, 234, 237],
    "DII":  [859, 862, 865, 868, 871, 874],
    "DIII": [1316, 1319, 1322, 1325, 1328],
    "DIV":  [1636, 1639, 1642, 1645, 1648, 1651],
}

# IFM inactivation motif (DIII–DIV linker). F1499 is the central
# residue. Receptor pocket residues are formed by C-terminal end of
# DIV S6 and the DIII-DIV linker; Phe1764 (DIV S6 helix end) is one
# of the canonical receptor residues.
SCN1A_IFM = ("F", 1499)
SCN1A_IFM_RECEPTOR = [("F", 1764), ("N", 1750), ("Y", 1762)]

# S6 bundle-crossing region (activation gate). Cytoplasmic ends of
# the four S6 helices form the activation gate. Approximate
# 12-residue windows from the AF3 model.
SCN1A_S6_GATE = {
    "DI":   range(420, 432),
    "DII":  range(975, 987),
    "DIII": range(1455, 1467),
    "DIV":  range(1769, 1781),
}

# Equivalent residues in cryo-EM Nav structures, using authoritative
# residue maps from Pan/Yan publications. Where chain assignments
# differ across PDBs we provide a single chain hint; the script
# tries the hint first, then falls back to chain auto-detection.
NAV_REFERENCE_MAPS: dict[str, dict] = {
    # Each PDB is annotated with:
    #   label            -- channel identity for figure legends
    #   published_state  -- conformational state from the publication
    #   alpha_chain      -- chain ID(s) of the α-subunit (others are
    #                       β-subunits, toxins, lipids); search DEKA
    #                       only in α-subunits
    # DEKA residues are auto-discovered by sequence motif search
    # within α-subunit chains.
    "7K48": {
        "label": "Nav1.7-NavAb chimera (resting)",
        "published_state": "resting",
        "alpha_chain": None,  # NavAb chimera, distinct numbering
        "skip_residue_specific": True,
    },
    "6J8E": {
        "label": "Nav1.2 + β2 + KIIIA (pore blocked)",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
    "7XVE": {
        "label": "Nav1.7 mutant class-I",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
    "6UZ3": {
        "label": "Nav1.5 cardiac",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
    "7DTC": {
        "label": "Nav1.5 E1784K",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
    "7XM9": {
        "label": "Nav1.7 + β1/β2 + XEN907",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
    "6AGF": {
        "label": "Nav1.4 + β1",
        "published_state": "fast-inactivated",
        "alpha_chain": ["A"],
    },
}


# Sequence motifs for the four DEKA selectivity-filter residues.
# Each pattern is conserved across mammalian Nav1.x isoforms.
# Group 1 of each regex captures the SF residue. The patterns
# below were verified against SCN1A (Nav1.1, P35498) and the
# mammalian Nav cryo-EM structures used as references.
DEKA_MOTIFS = {
    "DI":   (re.compile(r"T[QHN](D)FW"),                  "D"),  # T(QH)DFW
    "DII":  (re.compile(r"CG(E)WIE|G(E)WIE|MG(E)WIE"),    "E"),  # (C/M)GEWIE
    "DIII": (re.compile(r"[VAL][AS][TS]F(K)GW|F(K)GWMA"), "K"),  # ATFKGW
    "DIV":  (re.compile(r"T[ST](A)GW|TS(A)GW"),           "A"),  # TS(or T)AGW
}


# ----------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------
@dataclass
class Frame:
    """Local pore frame: origin = SF Cα centroid; +z = extracellular."""
    origin: np.ndarray   # 3-vector
    z_hat: np.ndarray    # unit vector along pore axis, pointing extracellular


def parse_structure(path: Path):
    parser = PDB.PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return parser.get_structure(path.stem, str(path))


def find_residue(structure, target_one_letter: str, resnum: int):
    """Return the first residue across all chains/models matching the
    given one-letter aa code and residue number, or None.

    AA matching uses one-letter; falls back to permissive match if
    three_to_index fails (some PDBs use HETATM for selenomet etc.)."""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # skip HETATM
                    continue
                if residue.id[1] != resnum:
                    continue
                resname = residue.get_resname()
                try:
                    one = index_to_one(three_to_index(resname))
                except (KeyError, IndexError):
                    one = None
                if one == target_one_letter:
                    return residue
        # Stop after first model
        break
    return None


def get_ca(residue) -> np.ndarray | None:
    if residue is None or "CA" not in residue:
        return None
    return np.asarray(residue["CA"].coord, dtype=float)


def build_frame(structure, deka_map: dict) -> Frame | None:
    """Build a pore-aligned local frame from DEKA Cα + S6 gate centroid."""
    deka_ca = []
    for _, (aa, num) in deka_map.items():
        ca = get_ca(find_residue(structure, aa, num))
        if ca is not None:
            deka_ca.append(ca)
    if len(deka_ca) < 3:
        return None
    sf_centroid = np.mean(deka_ca, axis=0)

    # Find any S6 gate residues to define +z direction (pore axis).
    # For SCN1A use SCN1A_S6_GATE; for cryo-EM PDBs caller must
    # provide gate map separately if numbering differs.
    return Frame(origin=sf_centroid, z_hat=np.array([0, 0, 1.0]))


def refine_z_axis(structure, frame: Frame, gate_map: dict) -> Frame:
    """Replace +z direction with the SF→gate vector (pore axis)."""
    gate_ca = []
    for _, residues in gate_map.items():
        for resnum in residues:
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                ca = get_ca(find_residue(structure, aa, resnum))
                if ca is not None:
                    gate_ca.append(ca)
                    break
    if len(gate_ca) < 4:
        return frame
    gate_centroid = np.mean(gate_ca, axis=0)
    axis = frame.origin - gate_centroid  # gate -> SF, so SF is +z (extracellular)
    norm = np.linalg.norm(axis)
    if norm < 1e-6:
        return frame
    return Frame(origin=frame.origin, z_hat=axis / norm)


def project_z(coord: np.ndarray, frame: Frame) -> float:
    return float(np.dot(coord - frame.origin, frame.z_hat))


def perpendicular_distance(coord: np.ndarray, frame: Frame) -> float:
    rel = coord - frame.origin
    along = np.dot(rel, frame.z_hat) * frame.z_hat
    perp = rel - along
    return float(np.linalg.norm(perp))


def pore_radius_profile(structure, frame: Frame,
                        z_min: float = -25.0, z_max: float = 25.0,
                        z_step: float = 1.0,
                        slab_half: float = 1.5,
                        chain_filter: list[str] | None = None) -> pd.DataFrame:
    """Minimum perpendicular distance from pore axis to any heavy atom
    in a slab around each z. Approximates the channel radius
    interpretable as the maximum sphere that fits at this z."""
    atoms = []
    for model in structure:
        for chain in model:
            if chain_filter and chain.id not in chain_filter:
                continue
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                for atom in residue:
                    if atom.element == "H":
                        continue
                    coord = np.asarray(atom.coord, dtype=np.float64)
                    if not np.all(np.isfinite(coord)):
                        continue
                    if np.linalg.norm(coord) > 1e4:  # discard absurd outliers
                        continue
                    atoms.append(coord)
        break
    if not atoms:
        return pd.DataFrame()
    coords = np.vstack(atoms).astype(np.float64)
    origin = np.asarray(frame.origin, dtype=np.float64)
    z_hat = np.asarray(frame.z_hat, dtype=np.float64)
    z_hat = z_hat / max(np.linalg.norm(z_hat), 1e-9)
    rel = coords - origin
    z = rel @ z_hat
    perp = np.linalg.norm(rel - np.outer(z, z_hat), axis=1)
    rows = []
    z_grid = np.arange(z_min, z_max + z_step, z_step)
    for zc in z_grid:
        mask = (z >= zc - slab_half) & (z <= zc + slab_half)
        rows.append({
            "z": float(zc),
            "min_radius_A": float(perp[mask].min()) if mask.any() else float("nan"),
            "n_atoms_in_slab": int(mask.sum()),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Per-structure analysis
# ----------------------------------------------------------------------
def analyze_scn1a() -> dict:
    print(f"[scn1a] parsing {AF3_PDB}", file=sys.stderr)
    structure = parse_structure(AF3_PDB)
    frame0 = build_frame(structure, SCN1A_DEKA)
    if frame0 is None:
        raise RuntimeError("SCN1A DEKA residues not found in AF3 model")
    frame = refine_z_axis(structure, frame0, SCN1A_S6_GATE)

    # DEKA Cα pairwise distances
    deka_ca = {dom: get_ca(find_residue(structure, aa, num))
               for dom, (aa, num) in SCN1A_DEKA.items()}
    deka_pairs = {}
    pair_keys = [("DI", "DIII"), ("DII", "DIV"), ("DI", "DII"),
                 ("DIII", "DIV"), ("DI", "DIV"), ("DII", "DIII")]
    for a, b in pair_keys:
        ca_a, ca_b = deka_ca[a], deka_ca[b]
        if ca_a is not None and ca_b is not None:
            deka_pairs[f"DEKA_{a}_{b}_CA"] = float(np.linalg.norm(ca_a - ca_b))

    # S4 outermost gating arginine z-position per VSD.
    # We project all S4 arginines onto the pore axis and report the
    # most extracellular (largest z) as R1 — this is robust against
    # uncertainty over the exact sequence position of R1.
    s4_z = {}
    for dom, positions in SCN1A_S4_ARGS.items():
        zs = []
        for num in positions:
            ca = get_ca(find_residue(structure, "R", num))
            if ca is not None:
                zs.append(project_z(ca, frame))
        s4_z[f"S4_R1_{dom}_z"] = max(zs) if zs else float("nan")
        s4_z[f"S4_zspan_{dom}"] = (max(zs) - min(zs)) if zs else float("nan")

    # IFM motif distance to receptor centroid
    ifm_ca = get_ca(find_residue(structure, *SCN1A_IFM))
    rec_cas = [get_ca(find_residue(structure, *r)) for r in SCN1A_IFM_RECEPTOR]
    rec_cas = [c for c in rec_cas if c is not None]
    if ifm_ca is not None and rec_cas:
        rec_centroid = np.mean(rec_cas, axis=0)
        ifm_dist = float(np.linalg.norm(ifm_ca - rec_centroid))
    else:
        ifm_dist = float("nan")

    # Pore radius profile
    profile = pore_radius_profile(structure, frame)
    profile["structure"] = "SCN1A_AF3"

    # Activation gate radius = min over z ∈ (-15, 0) Å (cytoplasmic side)
    gate_band = profile[(profile["z"] >= -15) & (profile["z"] <= 0)]
    gate_radius = float(gate_band["min_radius_A"].min()) if len(gate_band) else float("nan")

    return {
        "structure": "SCN1A_AF3",
        "label": "SCN1A AF3 (this study)",
        **deka_pairs,
        **s4_z,
        "IFM_to_receptor_A": ifm_dist,
        "activation_gate_radius_A": gate_radius,
        "_profile": profile,
    }


def chain_sequence(chain) -> tuple[str, list[int]]:
    """Return (single-letter sequence, ordered residue numbers) for a chain."""
    seq = []
    nums = []
    for r in chain:
        if r.id[0] != " ":
            continue
        try:
            one = index_to_one(three_to_index(r.get_resname()))
        except (KeyError, IndexError):
            one = "X"
        seq.append(one)
        nums.append(r.id[1])
    return "".join(seq), nums


def find_deka_in_chain(chain) -> dict[str, tuple[str, int]] | None:
    """Use DEKA_MOTIFS to locate the SF residue numbers in a chain.
    Returns {'DI': ('D', resnum), 'DII': ('E', resnum), ...} or None
    if any domain motif fails."""
    seq, nums = chain_sequence(chain)
    out = {}
    for dom, (regex, aa) in DEKA_MOTIFS.items():
        match = regex.search(seq)
        if match is None:
            return None
        # First non-empty group is the SF residue
        for i in range(1, regex.groups + 1):
            if match.group(i):
                idx = match.start(i)
                resnum = nums[idx]
                out[dom] = (aa, resnum)
                break
    return out


def analyze_reference(pdb_id: str, cfg: dict) -> dict:
    pdb_path = PDB_DIR / f"{pdb_id}.pdb"
    print(f"[ref ] parsing {pdb_path.name}: {cfg['label']}", file=sys.stderr)
    structure = parse_structure(pdb_path)
    out = {"structure": pdb_id, "label": cfg["label"],
           "published_state": cfg.get("published_state", "")}

    if cfg.get("skip_residue_specific"):
        # NavAb chimera or otherwise non-standard numbering: compute
        # pore profile aligned to the protein principal axis.
        cas = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " " or "CA" not in residue:
                        continue
                    cas.append(np.asarray(residue["CA"].coord, dtype=float))
            break
        if not cas:
            return out
        cas_arr = np.vstack(cas).astype(np.float64)
        origin = cas_arr.mean(axis=0)
        centered = cas_arr - origin
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        z_hat = vt[0].astype(np.float64)
        frame = Frame(origin=origin.astype(np.float64), z_hat=z_hat)
        profile = pore_radius_profile(structure, frame)
        profile["structure"] = pdb_id
        out["_profile"] = profile
        gate_band = profile[(profile["z"] >= -15) & (profile["z"] <= 0)]
        out["activation_gate_radius_A"] = (
            float(gate_band["min_radius_A"].min()) if len(gate_band) else float("nan")
        )
        return out

    # Auto-discover DEKA via sequence motifs. Search all candidate
    # α-subunit chains and pick the one yielding all four matches.
    deka_map = None
    found_chain = None
    candidate_chains = cfg.get("alpha_chain") or [c.id for m in structure for c in m for _ in [None]][:1]
    if candidate_chains is None:
        candidate_chains = [c.id for c in next(iter(structure))]
    for model in structure:
        for chain in model:
            if chain.id not in candidate_chains:
                continue
            deka_map = find_deka_in_chain(chain)
            if deka_map is not None:
                found_chain = chain.id
                break
        if deka_map is not None:
            break
    # Fallback: try every chain
    if deka_map is None:
        for model in structure:
            for chain in model:
                deka_map = find_deka_in_chain(chain)
                if deka_map is not None:
                    found_chain = chain.id
                    break
            if deka_map is not None:
                break
    if deka_map is None:
        print(f"[ref ]   WARN: DEKA motifs not found in {pdb_id}; "
              f"falling back to SVD-aligned pore profile", file=sys.stderr)
        # Pick the longest chain in the structure as the α-subunit
        # and align by its principal axis through its centroid.
        best_chain = None
        best_n = 0
        for model in structure:
            for chain in model:
                cas = [r for r in chain if r.id[0] == " " and "CA" in r]
                if len(cas) > best_n:
                    best_n = len(cas)
                    best_chain = chain
            break
        if best_chain is None:
            return out
        cas_arr = np.vstack([
            np.asarray(r["CA"].coord, dtype=np.float64) for r in best_chain
            if r.id[0] == " " and "CA" in r
        ])
        origin = cas_arr.mean(axis=0)
        centered = cas_arr - origin
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        z_hat = vt[0]
        frame = Frame(origin=origin, z_hat=z_hat)
        profile = pore_radius_profile(structure, frame, chain_filter=[best_chain.id])
        profile["structure"] = pdb_id
        out["_profile"] = profile
        gate_band = profile[(profile["z"] >= -25) & (profile["z"] <= -10)]
        out["activation_gate_radius_A"] = (
            float(gate_band["min_radius_A"].min()) if len(gate_band) else float("nan")
        )
        return out

    print(f"[ref ]   DEKA in chain {found_chain}: "
          f"D{deka_map['DI'][1]} E{deka_map['DII'][1]} K{deka_map['DIII'][1]} A{deka_map['DIV'][1]}",
          file=sys.stderr)

    # Restrict subsequent residue lookups to the discovered chain.
    target_chain = None
    for model in structure:
        for chain in model:
            if chain.id == found_chain:
                target_chain = chain
                break
        break

    def _ca_in_chain(aa: str, num: int) -> np.ndarray | None:
        if target_chain is None:
            return None
        for r in target_chain:
            if r.id[0] != " " or r.id[1] != num:
                continue
            try:
                one = index_to_one(three_to_index(r.get_resname()))
            except (KeyError, IndexError):
                one = None
            if one == aa and "CA" in r:
                return np.asarray(r["CA"].coord, dtype=float)
        return None

    # SF Cα centroid → frame origin
    deka_ca = {dom: _ca_in_chain(aa, num) for dom, (aa, num) in deka_map.items()}
    valid = [c for c in deka_ca.values() if c is not None]
    if len(valid) < 3:
        return out
    sf_centroid = np.mean(valid, axis=0)
    # Pore axis: SF→C-terminal centroid of the chain (cytoplasmic side).
    chain_cas = [np.asarray(r["CA"].coord, dtype=float)
                 for r in target_chain if r.id[0] == " " and "CA" in r]
    if len(chain_cas) < 100:
        return out
    chain_cas_arr = np.vstack(chain_cas)
    centered = chain_cas_arr - sf_centroid
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    z_hat = vt[0]
    # Orient so that +z points away from the cytoplasmic centroid
    # (extracellular). Use the bottom 25% of CAs in the SVD axis as
    # cytoplasmic anchor.
    proj = centered @ z_hat
    cyto_centroid = chain_cas_arr[proj < np.percentile(proj, 25)].mean(axis=0)
    if np.dot(sf_centroid - cyto_centroid, z_hat) < 0:
        z_hat = -z_hat
    frame = Frame(origin=sf_centroid, z_hat=z_hat)

    # DEKA Cα distances
    pair_keys = [("DI", "DIII"), ("DII", "DIV"), ("DI", "DII"),
                 ("DIII", "DIV"), ("DI", "DIV"), ("DII", "DIII")]
    for a, b in pair_keys:
        if deka_ca[a] is not None and deka_ca[b] is not None:
            out[f"DEKA_{a}_{b}_CA"] = float(np.linalg.norm(deka_ca[a] - deka_ca[b]))

    # Pore radius profile
    profile = pore_radius_profile(structure, frame, chain_filter=[found_chain])
    profile["structure"] = pdb_id
    out["_profile"] = profile
    gate_band = profile[(profile["z"] >= -15) & (profile["z"] <= -5)]
    out["activation_gate_radius_A"] = (
        float(gate_band["min_radius_A"].min()) if len(gate_band) else float("nan")
    )
    return out


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def make_figure(records: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plt.rcParams.update({"font.size": 10, "font.family": "Arial"})

    # (A) S4 R1 z-positions per VSD
    ax = axes[0, 0]
    rows = []
    for r in records:
        for dom in ("DI", "DII", "DIII", "DIV"):
            key = f"S4_R1_{dom}_z"
            if key in r and not np.isnan(r.get(key, float("nan"))):
                rows.append({"structure": r["structure"], "domain": dom,
                             "z": r[key], "label": r["label"]})
    df_s4 = pd.DataFrame(rows)
    if not df_s4.empty:
        pivot = df_s4.pivot(index="structure", columns="domain", values="z")
        pivot = pivot.reindex([r["structure"] for r in records if r["structure"] in pivot.index])
        pivot.plot(kind="bar", ax=ax, width=0.8,
                   color={"DI": "#1565C0", "DII": "#F57C00",
                          "DIII": "#388E3C", "DIV": "#D32F2F"})
        ax.set_ylabel("S4 outermost Arg Cα z\n(Å, +z = extracellular)", fontsize=10)
        ax.set_title("(A) S4 gating-charge z-position per VSD",
                     fontsize=11, fontweight="bold")
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.text(0.02, 0.97, "(activated)", transform=ax.transAxes, fontsize=8,
                color="#444444", va="top", style="italic")
        ax.text(0.02, 0.03, "(resting)", transform=ax.transAxes, fontsize=8,
                color="#444444", va="bottom", style="italic")
        ax.legend(title="VSD", fontsize=8, ncol=2, loc="upper right")
        ax.tick_params(axis="x", rotation=30, labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        ax.text(0.5, 0.5, "no S4 data", ha="center", va="center")

    # (B) DEKA Cα diagonal distances
    ax = axes[0, 1]
    diag1 = []
    diag2 = []
    labels = []
    for r in records:
        if "DEKA_DI_DIII_CA" in r and "DEKA_DII_DIV_CA" in r:
            labels.append(r["structure"])
            diag1.append(r["DEKA_DI_DIII_CA"])
            diag2.append(r["DEKA_DII_DIV_CA"])
    if labels:
        idx = np.arange(len(labels))
        w = 0.35
        ax.bar(idx - w / 2, diag1, w, label="DI–DIII (D–K)", color="#1565C0", alpha=0.8)
        ax.bar(idx + w / 2, diag2, w, label="DII–DIV (E–A)", color="#D32F2F", alpha=0.8)
        ax.set_xticks(idx)
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("Cα–Cα distance (Å)", fontsize=10)
        ax.set_title("(B) DEKA selectivity-filter Cα diagonals",
                     fontsize=11, fontweight="bold")
        ax.legend(fontsize=8, loc="lower right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        ax.text(0.5, 0.5, "no DEKA data", ha="center", va="center")

    # (C) Pore radius profile — emphasise SCN1A AF3 + key references
    ax = axes[1, 0]
    emphasis = {"SCN1A_AF3": ("#D32F2F", 2.5, 1.0, "SCN1A AF3 (this study)"),
                "7K48": ("#1565C0", 1.8, 0.95, "7K48 Nav1.7 resting"),
                "6J8E": ("#7B1FA2", 1.8, 0.95, "6J8E Nav1.2 inactivated"),
                "6UZ3": ("#388E3C", 1.4, 0.7, "6UZ3 Nav1.5"),
                "6AGF": ("#F57C00", 1.4, 0.7, "6AGF Nav1.4"),
                "7XVE": ("#90A4AE", 1.0, 0.5, "7XVE Nav1.7"),
                "7DTC": ("#90A4AE", 1.0, 0.5, "7DTC Nav1.5 mut"),
                "7XM9": ("#90A4AE", 1.0, 0.5, "7XM9 Nav1.7 + β1/2")}
    for r in records:
        prof = r.get("_profile")
        if prof is None or len(prof) == 0:
            continue
        styling = emphasis.get(r["structure"],
                               ("#BBBBBB", 0.8, 0.4, r["structure"]))
        color, lw, alpha, label = styling
        # 3-point rolling smooth on the radius profile to dampen noise
        rad = prof["min_radius_A"].rolling(3, center=True, min_periods=1).mean()
        ax.plot(rad, prof["z"], color=color, alpha=alpha, linewidth=lw, label=label)
    ax.set_xlabel("min radius from pore axis (Å)", fontsize=10)
    ax.set_ylabel("z (Å, 0 = SF mid-plane, +z = extracellular)", fontsize=10)
    ax.set_title("(C) Pore radius profile", fontsize=11, fontweight="bold")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--", alpha=0.6)
    # Shade the activation-gate band
    ax.axhspan(-15, 0, alpha=0.06, color="black")
    ax.text(13, -7.5, "activation\ngate band", ha="right", va="center", fontsize=7,
            color="#666666", style="italic")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.set_xlim(0, 15)
    ax.set_ylim(-25, 25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # (D) IFM-related panel
    ax = axes[1, 1]
    rows = []
    for r in records:
        if "IFM_to_receptor_A" in r and not math.isnan(r.get("IFM_to_receptor_A", float("nan"))):
            rows.append({"structure": r["structure"], "metric": "IFM-receptor (Å)",
                         "value": r["IFM_to_receptor_A"]})
        if "IFM_z_A" in r and not math.isnan(r.get("IFM_z_A", float("nan"))):
            rows.append({"structure": r["structure"], "metric": "IFM Cα z (Å)",
                         "value": r["IFM_z_A"]})
        if "activation_gate_radius_A" in r and not math.isnan(r.get("activation_gate_radius_A", float("nan"))):
            rows.append({"structure": r["structure"], "metric": "act-gate radius (Å)",
                         "value": r["activation_gate_radius_A"]})
    df_ifm = pd.DataFrame(rows)
    if not df_ifm.empty:
        # Show only activation-gate radius per structure (most robust comparator)
        gate_data = [(r["structure"], r.get("activation_gate_radius_A", float("nan")),
                      r.get("published_state", ""))
                     for r in records if not pd.isna(r.get("activation_gate_radius_A", float("nan")))]
        if gate_data:
            structures = [g[0] for g in gate_data]
            radii = [g[1] for g in gate_data]
            states = [g[2] for g in gate_data]
            colors = ["#D32F2F" if s == "SCN1A_AF3" else
                      "#1565C0" if "resting" in st else "#F57C00"
                      for s, st in zip(structures, states)]
            ax.bar(range(len(structures)), radii, color=colors, alpha=0.8)
            ax.set_xticks(range(len(structures)))
            ax.set_xticklabels(structures, rotation=30, ha="right", fontsize=8)
            ax.set_ylabel("Activation-gate min radius (Å)", fontsize=10)
            ax.axhline(2.5, color="#1565C0", linestyle="--", linewidth=0.8, alpha=0.6,
                       label="open-pore reference (~2.5 Å)")
            ax.set_title("(D) Activation-gate pore radius",
                         fontsize=11, fontweight="bold")
            ax.legend(fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.text(0.02, 0.97, "non-conducting if < 2 Å",
                    transform=ax.transAxes, fontsize=7.5,
                    color="#666666", va="top", style="italic")
    else:
        ax.text(0.5, 0.5, "no gate data", ha="center", va="center")

    fig.suptitle("Supplementary Figure S4. Conformational-state diagnostic of the SCN1A "
                 "AlphaFold3 model vs deposited cryo-EM Nav structures",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_base = FIG_DIR / "2605_Supplementary_Figure_S4_state_diagnostic"
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"[fig] wrote {out_base}.{{png,pdf,tiff}}", file=sys.stderr)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main() -> None:
    records = [analyze_scn1a()]
    for pdb_id, cfg in NAV_REFERENCE_MAPS.items():
        try:
            records.append(analyze_reference(pdb_id, cfg))
        except Exception as exc:  # noqa: BLE001
            print(f"[ref ]   ERROR {pdb_id}: {exc}", file=sys.stderr)

    # Wide-format diagnostic table
    drop_cols = {"_profile"}
    rows = [{k: v for k, v in r.items() if k not in drop_cols} for r in records]
    df = pd.DataFrame(rows).set_index("structure")
    df.to_csv(RESULTS_DIR / "structural_state_diagnostic.tsv", sep="\t")
    print(df.round(2).to_string(), file=sys.stderr)
    print(f"[out] {RESULTS_DIR / 'structural_state_diagnostic.tsv'}", file=sys.stderr)

    # Long-form pore profile
    profiles = [r["_profile"] for r in records if r.get("_profile") is not None]
    if profiles:
        full = pd.concat(profiles, ignore_index=True)
        full.to_csv(RESULTS_DIR / "pore_radius_profile.tsv", sep="\t", index=False)
        print(f"[out] {RESULTS_DIR / 'pore_radius_profile.tsv'}", file=sys.stderr)

    make_figure(records)


if __name__ == "__main__":
    main()
