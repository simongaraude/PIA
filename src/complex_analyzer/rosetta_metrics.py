"""Optional PyRosetta-based interface metrics."""

from __future__ import annotations

import logging
from typing import Dict

logger = logging.getLogger(__name__)

HAVE_PYROSETTA = False

try:
    import pyrosetta
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

    HAVE_PYROSETTA = True
except ImportError:
    pass


def init_pyrosetta() -> bool:
    """Initialise PyRosetta once. Returns True on success."""
    if not HAVE_PYROSETTA:
        return False
    try:
        pyrosetta.init("-mute all -ignore_unrecognized_res true", silent=True)
        return True
    except Exception as exc:
        logger.warning("PyRosetta init failed: %s", exc)
        return False


def compute_rosetta_metrics(
    pdb_path: str, binder_chain: str, receptor_chain: str
) -> Dict:
    """Run InterfaceAnalyzerMover and return a dict of Rosetta metrics."""
    if not HAVE_PYROSETTA:
        return {}
    try:
        pose = pyrosetta.pose_from_pdb(pdb_path)
        interface_str = f"{receptor_chain}_{binder_chain}"

        iam = InterfaceAnalyzerMover(interface_str)
        iam.set_pack_separated(True)
        iam.set_compute_interface_sc(True)
        iam.set_compute_interface_energy(True)
        iam.set_calc_dSASA(True)
        iam.set_calc_hbond_sasaE(True)
        iam.apply(pose)

        return {
            "I_sc": iam.get_interface_score(),
            "dG_separated": iam.get_separated_interface_energy(),
            "dSASA_int": iam.get_interface_delta_sasa(),
            "shape_complementarity": iam.get_interface_sc(),
            "n_hbonds_rosetta": iam.get_interface_delta_hbond_unsat(),
        }
    except Exception as exc:
        logger.warning("Rosetta analysis failed for %s: %s", pdb_path, exc)
        return {}


def compute_ddG_vs_parental(
    complex_pdb: str, parental_pdb: str, interface_str: str
) -> float:
    """ddG = dG(mutant) − dG(parental) via Rosetta."""
    if not HAVE_PYROSETTA:
        return float("nan")
    try:
        iam = InterfaceAnalyzerMover(interface_str)
        iam.set_pack_separated(True)

        pose_mut = pyrosetta.pose_from_pdb(complex_pdb)
        iam.apply(pose_mut)
        dG_mut = iam.get_separated_interface_energy()

        pose_par = pyrosetta.pose_from_pdb(parental_pdb)
        iam.apply(pose_par)
        dG_par = iam.get_separated_interface_energy()

        return dG_mut - dG_par
    except Exception as exc:
        logger.warning("ddG computation failed: %s", exc)
        return float("nan")
