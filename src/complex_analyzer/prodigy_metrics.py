"""PRODIGY-based binding affinity prediction.

Fully open-source (Apache 2.0), no registration or license required.
pip install prodigy-prot

PRODIGY predicts:
    - ΔG binding (kcal/mol)   — analogous to Rosetta dG_separated / FoldX interaction energy
    - Kd (M)                  — dissociation constant
    - Contact decomposition   — charged-charged, charged-polar, charged-apolar,
                                polar-polar, apolar-polar, apolar-apolar
    - NIS (non-interacting surface) — % apolar and % charged

Reference:
    Vangone & Bonvin, eLife 2015 (10.7554/eLife.07454)
    Xue et al., Bioinformatics 2016 (10.1093/bioinformatics/btw514)
"""

from __future__ import annotations

import io
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict

logger = logging.getLogger(__name__)

# Try importing prodigy as a Python library first
try:
    import prodigy_prot  # noqa: F401
    HAVE_PRODIGY_LIB = False  # Force CLI mode — lib API differs across versions
except ImportError:
    HAVE_PRODIGY_LIB = False

# Check if prodigy CLI is available
HAVE_PRODIGY_CLI = shutil.which("prodigy") is not None

HAVE_PRODIGY = HAVE_PRODIGY_LIB or HAVE_PRODIGY_CLI


def _ensure_pdb(pdb_path: str) -> str:
    """If input is CIF, convert to PDB for PRODIGY compatibility."""
    ext = os.path.splitext(pdb_path)[1].lower()
    if ext in (".cif", ".mmcif"):
        from Bio.PDB import PDBIO

        from .structure import parse_structure

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        struct, _model = parse_structure(pdb_path, "conv")
        io = PDBIO()
        io.set_structure(struct)
        io.save(tmp.name)
        return tmp.name
    return pdb_path


def compute_prodigy_metrics(
    pdb_path: str,
    chain_a: str = "A",
    chain_b: str = "B",
    temperature: float = 25.0,
) -> Dict[str, float]:
    """
    Run PRODIGY on a protein-protein complex.

    Returns dict with:
        - prodigy_dG          : predicted ΔG binding (kcal/mol), negative = favorable
        - prodigy_Kd          : predicted dissociation constant (M)
        - prodigy_n_contacts  : total intermolecular contacts
        - prodigy_cc          : charged-charged contacts
        - prodigy_cp          : charged-polar contacts
        - prodigy_ca          : charged-apolar contacts
        - prodigy_pp          : polar-polar contacts
        - prodigy_ap          : apolar-polar contacts
        - prodigy_aa          : apolar-apolar contacts
        - prodigy_nis_apolar  : % apolar NIS residues
        - prodigy_nis_charged : % charged NIS residues
    """
    if not HAVE_PRODIGY:
        return {}

    converted = _ensure_pdb(pdb_path)
    cleanup = converted != pdb_path

    try:
        if HAVE_PRODIGY_LIB:
            return _run_prodigy_lib(converted, chain_a, chain_b, temperature)
        else:
            return _run_prodigy_cli(converted, chain_a, chain_b, temperature)
    except Exception as exc:
        logger.warning("PRODIGY failed for %s: %s", pdb_path, exc)
        return {}
    finally:
        if cleanup and os.path.exists(converted):
            os.unlink(converted)


def _run_prodigy_cli(
    pdb_path: str,
    chain_a: str,
    chain_b: str,
    temperature: float,
) -> Dict[str, float]:
    """Run PRODIGY via command line and parse output."""
    cmd = [
        "prodigy",
        pdb_path,
        "--selection", chain_a, chain_b,
        "--temperature", str(temperature),
    ]

    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=120
    )

    output = result.stdout + result.stderr
    return _parse_prodigy_output(output)


def _run_prodigy_lib(
    pdb_path: str,
    chain_a: str,
    chain_b: str,
    temperature: float,
) -> Dict[str, float]:
    """Run PRODIGY via Python library (prodigy_prot)."""
    old_stdout = sys.stdout
    sys.stdout = buffer = io.StringIO()

    try:
        from prodigy_prot.modules.prodigy import Prodigy
        structure = Prodigy(pdb_path, selection=[[chain_a], [chain_b]],
                           temperature=temperature)
        structure.predict()
        structure.print_prediction()
        output = buffer.getvalue()
    finally:
        sys.stdout = old_stdout

    return _parse_prodigy_output(output)


def _parse_prodigy_output(output: str) -> Dict[str, float]:
    """Parse PRODIGY text output into a metrics dict."""
    metrics: Dict[str, float] = {}

    for line in output.splitlines():
        line = line.strip()

        if "intermolecular contacts:" in line.lower():
            try:
                metrics["prodigy_n_contacts"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "charged-charged" in line.lower():
            try:
                metrics["prodigy_cc"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "charged-polar" in line.lower():
            try:
                metrics["prodigy_cp"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "charged-apolar" in line.lower():
            try:
                metrics["prodigy_ca"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "polar-polar" in line.lower():
            try:
                metrics["prodigy_pp"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "apolar-polar" in line.lower():
            try:
                metrics["prodigy_ap"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "apolar-apolar" in line.lower():
            try:
                metrics["prodigy_aa"] = int(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "apolar nis" in line.lower():
            try:
                metrics["prodigy_nis_apolar"] = float(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "charged nis" in line.lower():
            try:
                metrics["prodigy_nis_charged"] = float(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "binding affinity" in line.lower():
            try:
                val = line.split(":")[-1].strip()
                # Remove units like "kcal.mol-1"
                metrics["prodigy_dG"] = float(val.split()[0])
            except (ValueError, IndexError):
                pass
        elif "dissociation constant" in line.lower():
            try:
                val = line.split(":")[-1].strip()
                # Parse scientific notation like "1.3e-07"
                metrics["prodigy_Kd"] = float(val.split()[0])
            except (ValueError, IndexError):
                pass

    return metrics


def compute_prodigy_ddG(
    complex_pdb: str,
    parental_pdb: str,
    chain_a: str = "A",
    chain_b: str = "B",
    temperature: float = 25.0,
) -> float:
    """ΔΔG = ΔG(design) - ΔG(parental). Negative = better binding than parental."""
    m_design = compute_prodigy_metrics(complex_pdb, chain_a, chain_b, temperature)
    m_parent = compute_prodigy_metrics(parental_pdb, chain_a, chain_b, temperature)

    dg_d = m_design.get("prodigy_dG", float("nan"))
    dg_p = m_parent.get("prodigy_dG", float("nan"))

    try:
        return round(dg_d - dg_p, 4)
    except (TypeError, ValueError):
        return float("nan")
