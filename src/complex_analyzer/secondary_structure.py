"""Secondary structure assignment and fold classification."""

from __future__ import annotations

import logging
import os
import tempfile
from typing import Dict

import numpy as np
from Bio.PDB import DSSP, PDBIO

from .constants import HYDROPHOBIC_AA
from .structure import get_chain_residues

logger = logging.getLogger(__name__)

HELIX_CODES = {"H", "G", "I"}
SHEET_CODES = {"E", "B"}


def _ensure_pdb_for_dssp(pdb_path: str, model) -> str:
    """If input is CIF, write a temporary PDB for DSSP. Returns path to use."""
    ext = os.path.splitext(pdb_path)[1].lower()
    if ext in (".cif", ".mmcif"):
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        io = PDBIO()
        io.set_structure(model.get_parent())
        io.save(tmp.name)
        return tmp.name
    return pdb_path


def compute_secondary_structure(
    pdb_path: str, model, chain_id: str
) -> Dict[str, int]:
    """Assign SS via DSSP (falls back to CA heuristic)."""
    dssp_path = _ensure_pdb_for_dssp(pdb_path, model)
    cleanup = dssp_path != pdb_path

    try:
        try:
            dssp = DSSP(model, dssp_path, dssp="mkdssp")
        except Exception:
            dssp = DSSP(model, dssp_path, dssp="dssp")

        result = {"helix": 0, "sheet": 0, "loop": 0, "total": 0}
        for key in dssp.keys():
            if key[0] != chain_id:
                continue
            ss = dssp[key][2]
            result["total"] += 1
            if ss in HELIX_CODES:
                result["helix"] += 1
            elif ss in SHEET_CODES:
                result["sheet"] += 1
            else:
                result["loop"] += 1
        return result
    except Exception as exc:
        logger.warning("DSSP unavailable (%s); using CA heuristic.", exc)
        return _estimate_ss_from_ca(model, chain_id)
    finally:
        if cleanup and os.path.exists(dssp_path):
            os.unlink(dssp_path)


def _estimate_ss_from_ca(model, chain_id: str) -> Dict[str, int]:
    """Rough SS estimate from CA–CA distances when DSSP is absent."""
    residues = get_chain_residues(model, chain_id)
    ca = [
        r["CA"].get_vector().get_array() if "CA" in r else None
        for r in residues
    ]
    result = {"helix": 0, "sheet": 0, "loop": 0, "total": len(residues)}
    for i, coord in enumerate(ca):
        if coord is None:
            result["loop"] += 1
            continue
        if i + 3 < len(ca) and ca[i + 3] is not None:
            d = np.linalg.norm(np.array(coord) - np.array(ca[i + 3]))
            if 4.8 < d < 6.0:
                result["helix"] += 1
                continue
        result["loop"] += 1
    return result


def detect_folds(ss: Dict[str, int], sequence: str) -> str:
    """
    Simple fold classification from SS content.

    For proper fold assignment use Foldseek against CATH/SCOPe.
    """
    total = max(ss["total"], 1)
    pct_h = ss["helix"] / total
    pct_s = ss["sheet"] / total

    folds: list[str] = []
    if pct_h > 0.6:
        folds.append("all-alpha")
    elif pct_s > 0.4:
        folds.append("all-beta")
    elif pct_h > 0.25 and pct_s > 0.15:
        folds.append("alpha/beta")
    elif pct_h > 0.15 and pct_s > 0.1:
        folds.append("alpha+beta")
    else:
        folds.append("irregular/coil-rich")

    # Coiled-coil heuristic
    if pct_h > 0.5 and len(sequence) > 20:
        hydro = [1 if c in HYDROPHOBIC_AA else 0 for c in sequence]
        matches = sum(
            1
            for i in range(len(hydro) - 7)
            if hydro[i] and (hydro[i + 3] or hydro[i + 4])
        )
        if matches / max(len(sequence) - 7, 1) > 0.4:
            folds.append("possible-coiled-coil")

    # Disulfide-rich mini-protein
    if sequence.count("C") / max(len(sequence), 1) > 0.06 and len(sequence) < 80:
        folds.append("disulfide-rich-mini-protein")

    return "; ".join(folds)
