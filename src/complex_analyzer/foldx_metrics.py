"""FoldX-based interface energy metrics.

Requires:
    - FoldX binary (free for academics: https://foldxsuite.crg.eu)
    - rotabase.txt in the FoldX directory

FoldX AnalyseComplex computes:
    - Interaction Energy (ΔG_binding) — analogous to Rosetta I_sc / dG_separated
    - Solvation Hydrophobic + Solvation Polar — desolvation upon binding
    - Individual energy term decomposition (vdW, H-bonds, electrostatics, etc.)

For ddG vs parental: run AnalyseComplex on both, take the difference.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Default location on the cluster
FOLDX_DEFAULT_PATHS = [
    "/shared/foldx/foldx",
    "/shared/foldx/FoldX",
    "/usr/local/bin/foldx",
    os.path.expanduser("~/foldx/foldx"),
]


def find_foldx_binary(custom_path: Optional[str] = None) -> Optional[str]:
    """Locate the FoldX binary."""
    if custom_path and os.path.isfile(custom_path):
        return custom_path
    for p in FOLDX_DEFAULT_PATHS:
        if os.path.isfile(p):
            return p
    # Try PATH
    result = shutil.which("foldx")
    if result:
        return result
    return None


def find_rotabase(foldx_bin: str) -> Optional[str]:
    """Find rotabase.txt near the FoldX binary."""
    foldx_dir = os.path.dirname(foldx_bin)
    for candidate in [
        os.path.join(foldx_dir, "rotabase.txt"),
        os.path.join(foldx_dir, "..", "rotabase.txt"),
        os.path.join(os.getcwd(), "rotabase.txt"),
    ]:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def cif_to_pdb(cif_path: str, output_dir: str) -> str:
    """Convert CIF to PDB using BioPython. Returns path to PDB file."""
    from Bio.PDB import PDBIO

    from .structure import parse_structure

    basename = Path(cif_path).stem + ".pdb"
    pdb_path = os.path.join(output_dir, basename)

    struct, model = parse_structure(cif_path, "conv")
    io = PDBIO()
    io.set_structure(struct)
    io.save(pdb_path)
    return pdb_path


def run_repair_pdb(foldx_bin: str, pdb_path: str, work_dir: str) -> str:
    """Run FoldX RepairPDB. Returns path to repaired PDB."""
    pdb_name = os.path.basename(pdb_path)
    stem = Path(pdb_name).stem

    cmd = [
        foldx_bin,
        "--command=RepairPdb",
        f"--pdb={pdb_name}",
        f"--output-dir={work_dir}",
    ]

    result = subprocess.run(
        cmd, cwd=work_dir, capture_output=True, text=True, timeout=300
    )

    repaired = os.path.join(work_dir, f"{stem}_Repair.pdb")
    if not os.path.isfile(repaired):
        # Some FoldX versions use different naming
        for f in os.listdir(work_dir):
            if "Repair" in f and f.endswith(".pdb"):
                repaired = os.path.join(work_dir, f)
                break

    if not os.path.isfile(repaired):
        logger.warning("RepairPDB failed for %s: %s", pdb_path, result.stderr[:200])
        return pdb_path  # Fall back to unrepaired

    return repaired


def run_analyse_complex(
    foldx_bin: str,
    pdb_path: str,
    chain_a: str,
    chain_b: str,
    work_dir: str,
) -> Dict[str, float]:
    """
    Run FoldX AnalyseComplex and parse results.

    Returns dict with:
        - interaction_energy (ΔG_binding, kcal/mol)
        - backbone_hbond, sidechain_hbond, van_der_waals, electrostatics
        - solvation_polar, solvation_hydrophobic
        - entropy_sidechain, entropy_mainchain
        - total_stability (from Summary)
    """
    pdb_name = os.path.basename(pdb_path)

    cmd = [
        foldx_bin,
        "--command=AnalyseComplex",
        f"--pdb={pdb_name}",
        f"--analyseComplexChains={chain_a},{chain_b}",
        f"--output-dir={work_dir}",
    ]

    result = subprocess.run(
        cmd, cwd=work_dir, capture_output=True, text=True, timeout=300
    )

    metrics: Dict[str, float] = {}

    # Parse Interaction_*.fxout
    interaction_file = None
    for f in os.listdir(work_dir):
        if f.startswith("Interaction_") and f.endswith(".fxout"):
            interaction_file = os.path.join(work_dir, f)
            break

    if interaction_file and os.path.isfile(interaction_file):
        metrics.update(_parse_foldx_energy_file(interaction_file))

    # Parse Summary_*.fxout
    summary_file = None
    for f in os.listdir(work_dir):
        if f.startswith("Summary_") and f.endswith(".fxout"):
            summary_file = os.path.join(work_dir, f)
            break

    if summary_file and os.path.isfile(summary_file):
        summary = _parse_foldx_summary(summary_file)
        metrics.update(summary)

    if not metrics:
        logger.warning(
            "FoldX AnalyseComplex produced no output for %s: %s",
            pdb_path,
            result.stderr[:200],
        )

    return metrics


def _parse_foldx_energy_file(filepath: str) -> Dict[str, float]:
    """
    Parse FoldX Interaction_*.fxout or Indiv_energies_*.fxout.

    These are tab-separated with a header row starting with 'Pdb'.
    Columns (FoldX 5.x):
        Pdb, Group1, Group2, IntraclashesGroup1, IntraclashesGroup2,
        Backbone Hbond, Sidechain Hbond, Van der Waals, Electrostatics,
        Solvation Polar, Solvation Hydrophobic, Van der Waals clashes,
        entropy sidechain, entropy mainchain, sloop_entropy,
        mloop_entropy, cis_bond, torsional clash, backbone clash,
        helix dipole, water bridge, disulfide, electrostatic kon,
        partial covalent bonds, energy ionisation, Entropy Complex,
        Interface Energy (= Interaction Energy)
    """
    result: Dict[str, float] = {}

    try:
        with open(filepath) as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("///")]

        # Find header line
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith("Pdb"):
                header_idx = i
                break

        if header_idx is None:
            return result

        headers = lines[header_idx].split("\t")
        # Take the last data line (there may be multiple)
        data_line = lines[-1].split("\t")

        if len(data_line) < len(headers):
            return result

        col_map = {
            "Backbone Hbond": "foldx_backbone_hbond",
            "Sidechain Hbond": "foldx_sidechain_hbond",
            "Van der Waals": "foldx_van_der_waals",
            "Electrostatics": "foldx_electrostatics",
            "Solvation Polar": "foldx_solvation_polar",
            "Solvation Hydrophobic": "foldx_solvation_hydrophobic",
            "Van der Waals clashes": "foldx_vdw_clashes",
            "entropy sidechain": "foldx_entropy_sidechain",
            "entropy mainchain": "foldx_entropy_mainchain",
            "water bridge": "foldx_water_bridge",
            "disulfide": "foldx_disulfide",
            "electrostatic kon": "foldx_electrostatic_kon",
            "Entropy Complex": "foldx_entropy_complex",
            "IntraclashesGroup1": "foldx_intraclashes_receptor",
            "IntraclashesGroup2": "foldx_intraclashes_binder",
        }

        for i, header in enumerate(headers):
            header_clean = header.strip()
            if header_clean in col_map and i < len(data_line):
                try:
                    result[col_map[header_clean]] = float(data_line[i])
                except (ValueError, IndexError):
                    pass

        # Compute interaction energy (sum of all interface terms)
        # The last numeric column is typically the total interaction energy
        # but let's also compute it explicitly
        energy_terms = [
            "foldx_backbone_hbond", "foldx_sidechain_hbond",
            "foldx_van_der_waals", "foldx_electrostatics",
            "foldx_solvation_polar", "foldx_solvation_hydrophobic",
            "foldx_vdw_clashes", "foldx_entropy_sidechain",
            "foldx_entropy_mainchain", "foldx_water_bridge",
            "foldx_electrostatic_kon", "foldx_entropy_complex",
        ]
        total = sum(result.get(t, 0) for t in energy_terms)
        result["foldx_interaction_energy"] = round(total, 4)

        # dSASA proxy: solvation terms reflect desolvation upon binding
        desolv = result.get("foldx_solvation_hydrophobic", 0) + \
                 result.get("foldx_solvation_polar", 0)
        result["foldx_desolvation"] = round(desolv, 4)

    except Exception as exc:
        logger.warning("Failed to parse %s: %s", filepath, exc)

    return result


def _parse_foldx_summary(filepath: str) -> Dict[str, float]:
    """Parse Summary_*.fxout for complex stability."""
    result: Dict[str, float] = {}
    try:
        with open(filepath) as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("///")]

        for line in lines:
            if line.startswith("Pdb"):
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                try:
                    result["foldx_complex_stability"] = float(parts[-1])
                except (ValueError, IndexError):
                    pass
    except Exception as exc:
        logger.warning("Failed to parse summary %s: %s", filepath, exc)

    return result


def analyse_complex_foldx(
    pdb_path: str,
    chain_a: str = "A",
    chain_b: str = "B",
    foldx_bin: Optional[str] = None,
    repair: bool = True,
) -> Dict[str, float]:
    """
    Full FoldX analysis pipeline for one complex.

    1. Convert CIF→PDB if needed
    2. RepairPDB (optional but recommended)
    3. AnalyseComplex
    4. Return parsed metrics
    """
    foldx = find_foldx_binary(foldx_bin)
    if foldx is None:
        return {}

    rotabase = find_rotabase(foldx)

    work_dir = tempfile.mkdtemp(prefix="foldx_")
    try:
        # Convert CIF if needed
        ext = os.path.splitext(pdb_path)[1].lower()
        if ext in (".cif", ".mmcif"):
            pdb_file = cif_to_pdb(pdb_path, work_dir)
        else:
            pdb_file = os.path.join(work_dir, os.path.basename(pdb_path))
            shutil.copy2(pdb_path, pdb_file)

        # Copy rotabase.txt if found
        if rotabase:
            shutil.copy2(rotabase, os.path.join(work_dir, "rotabase.txt"))

        # Repair
        if repair:
            pdb_file = run_repair_pdb(foldx, pdb_file, work_dir)

        # AnalyseComplex
        metrics = run_analyse_complex(foldx, pdb_file, chain_a, chain_b, work_dir)

        return metrics

    except Exception as exc:
        logger.warning("FoldX pipeline failed for %s: %s", pdb_path, exc)
        return {}
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def compute_foldx_ddG(
    complex_pdb: str,
    parental_pdb: str,
    chain_a: str = "A",
    chain_b: str = "B",
    foldx_bin: Optional[str] = None,
) -> float:
    """ddG = interaction_energy(mutant) - interaction_energy(parental)."""
    metrics_mut = analyse_complex_foldx(complex_pdb, chain_a, chain_b, foldx_bin)
    metrics_par = analyse_complex_foldx(parental_pdb, chain_a, chain_b, foldx_bin)

    ie_mut = metrics_mut.get("foldx_interaction_energy", float("nan"))
    ie_par = metrics_par.get("foldx_interaction_energy", float("nan"))

    try:
        return round(ie_mut - ie_par, 4)
    except (TypeError, ValueError):
        return float("nan")
