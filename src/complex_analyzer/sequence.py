"""Sequence-level metrics: composition, charge, amphipathicity, cell penetrance."""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, Tuple

import numpy as np

from .constants import (
    AA_CHARGE,
    EISENBERG_HYDROPHOBICITY,
    HYDROPHOBIC_AA,
    KD_HYDROPHOBICITY,
    NEG_CHARGED_AA,
    POLAR_AA,
    POS_CHARGED_AA,
    STANDARD_AA,
)


def compute_composition(sequence: str) -> Dict:
    """Counts and percentages for charge classes and individual amino acids."""
    n = max(len(sequence), 1)
    aa_counts = Counter(sequence)

    n_pos = sum(aa_counts.get(a, 0) for a in POS_CHARGED_AA)
    n_neg = sum(aa_counts.get(a, 0) for a in NEG_CHARGED_AA)
    n_charged = n_pos + n_neg
    n_hydro = sum(aa_counts.get(a, 0) for a in HYDROPHOBIC_AA)
    n_polar = sum(aa_counts.get(a, 0) for a in POLAR_AA)

    net_charge = round(
        sum(AA_CHARGE.get(a, 0) * aa_counts.get(a, 0) for a in AA_CHARGE), 1
    )

    result = {
        "n_residues": len(sequence),
        "net_charge": net_charge,
        "n_charged": n_charged,
        "pct_charged": round(100 * n_charged / n, 1),
        "n_pos_charged": n_pos,
        "pct_pos_charged": round(100 * n_pos / n, 1),
        "n_neg_charged": n_neg,
        "pct_neg_charged": round(100 * n_neg / n, 1),
        "n_hydrophobic": n_hydro,
        "pct_hydrophobic": round(100 * n_hydro / n, 1),
        "n_polar": n_polar,
        "pct_polar": round(100 * n_polar / n, 1),
    }

    for aa in STANDARD_AA:
        result[f"count_{aa}"] = aa_counts.get(aa, 0)
        result[f"pct_{aa}"] = round(100 * aa_counts.get(aa, 0) / n, 1)

    return result


def compute_hydrophobic_moment(
    sequence: str, window: int = 11, angle: float = 100.0
) -> Tuple[float, float]:
    """
    Mean and max hydrophobic moment (Eisenberg scale).

    *angle* = 100° for α-helix, 160° for β-strand.
    """
    if not sequence:
        return 0.0, 0.0

    angle_rad = math.radians(angle)

    def _moment(subseq: str) -> float:
        h = [EISENBERG_HYDROPHOBICITY.get(a, 0) for a in subseq]
        n = len(h)
        sin_sum = sum(v * math.sin(i * angle_rad) for i, v in enumerate(h))
        cos_sum = sum(v * math.cos(i * angle_rad) for i, v in enumerate(h))
        return math.sqrt(sin_sum ** 2 + cos_sum ** 2) / n

    if len(sequence) < window:
        m = _moment(sequence)
        return m, m

    moments = [_moment(sequence[i : i + window]) for i in range(len(sequence) - window + 1)]
    return float(np.mean(moments)), float(np.max(moments))


def predict_cell_penetrance(
    sequence: str,
    net_charge: float,
    mean_hydrophobic_moment: float,
    pct_hydrophobic: float,
) -> float:
    """
    Heuristic cell-penetrance score (0–1).

    Combines features correlated with membrane translocation (inspired by
    CPPpred / CellPPD).  **Not a validated predictor** — for production use
    submit sequences to CPPpred, CellPPD, or MLCPP.
    """
    n = max(len(sequence), 1)
    score = 0.0

    # Positive charge (CPPs typically +4 to +12)
    if net_charge > 0:
        score += 0.25 * min(net_charge / 8.0, 1.0)

    # Arginine fraction
    score += 0.20 * min(sequence.count("R") / (0.25 * n), 1.0)

    # Amphipathicity
    score += 0.20 * min(mean_hydrophobic_moment / 0.5, 1.0)

    # Moderate hydrophobicity
    mean_kd = np.mean([KD_HYDROPHOBICITY.get(a, 0) for a in sequence]) if sequence else 0
    score += 0.15 * max(0.0, 1.0 - abs(mean_kd - 0.5) / 3.0)

    # Length penalty
    if n <= 30:
        score += 0.10
    elif n <= 60:
        score += 0.05
    elif n <= 100:
        score += 0.02
    else:
        score += 0.005

    # Aromatic anchoring (Trp/Phe)
    aromatic = (sequence.count("W") + sequence.count("F")) / n
    score += 0.10 * min(aromatic / 0.15, 1.0)

    return round(min(max(score, 0.0), 1.0), 4)
