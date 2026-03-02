"""Command-line interface for complex-analyzer."""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
import traceback

import numpy as np
import pandas as pd

from .analyzer import analyze_complex
from .constants import STANDARD_AA
from .rosetta_metrics import HAVE_PYROSETTA, init_pyrosetta

logger = logging.getLogger(__name__)


# Column ordering for the output CSV
_PRIORITY_COLS = [
    "filename", "binder_chain", "receptor_chain", "n_residues",
    # Docking
    "I_sc", "dI_sc", "dG_separated", "dSASA_int", "desolvation_energy",
    # Interface
    "buried_surface_area", "shape_complementarity", "n_interface_contacts",
    "interface_packing_density", "n_hbonds", "n_salt_bridges",
    "hydrophobic_interface_area",
    "n_interface_residues_binder", "n_interface_residues_receptor",
    # Hotspots
    "hotspot_residues_binder",
    # SASA
    "total_sasa", "binder_sasa", "receptor_sasa",
    # Charge & composition
    "net_charge", "n_charged", "pct_charged",
    "n_pos_charged", "pct_pos_charged", "n_neg_charged", "pct_neg_charged",
    "n_hydrophobic", "pct_hydrophobic", "n_polar", "pct_polar",
    # SS
    "n_helix_residues", "pct_helix", "n_sheet_residues", "pct_sheet",
    "n_loop_residues", "pct_loop",
    # Folds / misc
    "detected_folds", "predicted_cell_penetrance_score",
    "mean_hydrophobic_moment", "max_hydrophobic_moment",
    # FoldX decomposed terms
    "foldx_interaction_energy", "foldx_desolvation",
    "foldx_backbone_hbond", "foldx_sidechain_hbond",
    "foldx_van_der_waals", "foldx_electrostatics",
    "foldx_solvation_polar", "foldx_solvation_hydrophobic",
    "foldx_vdw_clashes", "foldx_entropy_sidechain", "foldx_entropy_mainchain",
    "foldx_water_bridge", "foldx_electrostatic_kon", "foldx_entropy_complex",
    "foldx_complex_stability",
    "foldx_intraclashes_receptor", "foldx_intraclashes_binder",
    # PRODIGY binding affinity
    "prodigy_dG", "prodigy_n_contacts",
    "prodigy_cc", "prodigy_cp", "prodigy_ca",
    "prodigy_pp", "prodigy_ap", "prodigy_aa",
    "prodigy_nis_apolar", "prodigy_nis_charged",
]


# Columns excluded from output
_EXCLUDE_COLS = {"filter_rmsd", "ddG_vs_parental", "prodigy_Kd"}


def _ordered_columns(df: pd.DataFrame) -> list[str]:
    cols = list(_PRIORITY_COLS)
    for aa in STANDARD_AA:
        cols.extend([f"count_{aa}", f"pct_{aa}"])
    ordered = [c for c in cols if c in df.columns and c not in _EXCLUDE_COLS]
    remaining = [c for c in df.columns if c not in ordered and c not in _EXCLUDE_COLS]
    return ordered + remaining


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    ap = argparse.ArgumentParser(
        description="Comprehensive structural analysis of protein-protein complexes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples
--------
  # Analyse all PDBs in ~/Desktop/all_structures (default)
  complex-analyzer

  # With parental and endogenous references
  complex-analyzer --parental_pdb parental.pdb --endogenous_pdb endogenous.pdb

  # Custom input directory
  complex-analyzer --input_dir ./other_complexes/ -o results.csv

  # Single file
  complex-analyzer --input_files complex1.pdb -o results.csv
""",
    )
    default_input = os.path.join(os.path.expanduser("~"), "Desktop", "all_structures")
    # On the Slurm cluster, structures live in /shared/all_structures
    if os.path.isdir("/shared/all_structures"):
        default_input = "/shared/all_structures"
    default_output = os.path.join(default_input, "complex_metrics.csv")
    ap.add_argument(
        "--input_dir",
        default=default_input,
        help="Directory of PDB files (default: ~/Desktop/all_structures)",
    )
    ap.add_argument("--input_files", nargs="+", help="Specific PDB file(s)")
    ap.add_argument(
        "-o", "--output",
        default=default_output,
        help="Output CSV (default: ~/Desktop/all_structures/complex_metrics.csv)",
    )
    ap.add_argument("-b", "--binder_chain", default="B")
    ap.add_argument("-r", "--receptor_chain", default="A")
    ap.add_argument("--parental_pdb", default=None)
    ap.add_argument("--endogenous_pdb", default=None)
    ap.add_argument("--endogenous_binder_chain", default=None)
    ap.add_argument("--endogenous_receptor_chain", default=None)
    ap.add_argument(
        "--foldx_bin", default=None,
        help="Path to FoldX binary (auto-detected from /shared/foldx/ or PATH)",
    )
    ap.add_argument(
        "--foldx_repair", action="store_true",
        help="Run FoldX RepairPDB before AnalyseComplex (slower but more accurate)",
    )
    args = ap.parse_args(argv)

    # Init Rosetta
    if HAVE_PYROSETTA:
        init_pyrosetta()
        logger.info("PyRosetta initialised.")
    else:
        logger.warning(
            "PyRosetta not found — Rosetta metrics (I_sc, dG, ddG, Sc) will be N/A."
        )

    # Check FoldX
    from .foldx_metrics import find_foldx_binary
    foldx_path = find_foldx_binary(args.foldx_bin)
    if foldx_path:
        logger.info("FoldX found at %s — will compute interaction energies.", foldx_path)
    elif not HAVE_PYROSETTA:
        logger.warning(
            "Neither PyRosetta nor FoldX found — I_sc, dG, dSASA, ddG will be N/A.\n"
            "Install FoldX from https://foldxsuite.crg.eu (free for academics)\n"
            "and place binary at /shared/foldx/foldx or pass --foldx_bin."
        )

    # Check PRODIGY (open-source binding affinity)
    from .prodigy_metrics import HAVE_PRODIGY
    if HAVE_PRODIGY:
        logger.info("PRODIGY found — will predict binding affinity (ΔG, Kd).")
    elif not HAVE_PYROSETTA and not foldx_path:
        logger.warning(
            "PRODIGY not found — install with: pip install prodigy-prot"
        )

    # Gather files
    pdb_files: list[str] = []
    if args.input_dir:
        if not os.path.isdir(args.input_dir):
            logger.error(
                "Input directory does not exist: %s\n"
                "Please create ~/Desktop/all_structures/ and place your PDB files there,\n"
                "or specify a different directory with --input_dir.",
                args.input_dir,
            )
            sys.exit(1)
        for ext in ("*.pdb", "*.cif", "*.mmcif", "*.ent"):
            pdb_files += sorted(glob.glob(os.path.join(args.input_dir, ext)))
    if args.input_files:
        pdb_files += args.input_files

    if not pdb_files:
        logger.error("No input files. Use --input_dir or --input_files.")
        sys.exit(1)

    logger.info("Analysing %d complex(es) …", len(pdb_files))

    rows: list[dict] = []
    for i, path in enumerate(pdb_files, 1):
        logger.info("[%d/%d] %s", i, len(pdb_files), os.path.basename(path))
        try:
            rows.append(
                analyze_complex(
                    pdb_path=path,
                    binder_chain=args.binder_chain,
                    receptor_chain=args.receptor_chain,
                    parental_pdb=args.parental_pdb,
                    endogenous_pdb=args.endogenous_pdb,
                    endogenous_binder_chain=args.endogenous_binder_chain,
                    endogenous_receptor_chain=args.endogenous_receptor_chain,
                    foldx_bin=args.foldx_bin,
                    foldx_repair=args.foldx_repair,
                )
            )
        except Exception as exc:
            logger.error("Failed %s: %s", path, exc)
            traceback.print_exc()
            rows.append({"filename": os.path.basename(path), "error": str(exc)})

    if not rows:
        logger.error("No results.")
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df = pd.DataFrame(rows)
    df = df[_ordered_columns(df)]
    df.to_csv(args.output, index=False, float_format="%.4f")
    logger.info("Written %d complexes × %d columns → %s", len(df), len(df.columns), args.output)

    # Summary
    summary_cols = [
        c for c in [
            "I_sc", "dG_separated", "buried_surface_area",
            "shape_complementarity", "n_hbonds", "n_salt_bridges",
            "predicted_cell_penetrance_score", "net_charge",
            "mean_hydrophobic_moment",
        ]
        if c in df.select_dtypes(include=[np.number]).columns
    ]
    if summary_cols:
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print(df[summary_cols].describe().round(3).to_string())
        print("=" * 70)
