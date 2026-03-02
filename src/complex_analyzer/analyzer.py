"""Core analysis pipeline — orchestrates all metric computations."""

from __future__ import annotations

import logging
import math
import os
from typing import Dict, Optional

from .foldx_metrics import (
    analyse_complex_foldx,
    compute_foldx_ddG,
    find_foldx_binary,
)
from .prodigy_metrics import (
    HAVE_PRODIGY,
    compute_prodigy_ddG,
    compute_prodigy_metrics,
)
from .rosetta_metrics import (
    HAVE_PYROSETTA,
    compute_ddG_vs_parental,
    compute_rosetta_metrics,
)
from .secondary_structure import compute_secondary_structure, detect_folds
from .sequence import compute_composition, compute_hydrophobic_moment, predict_cell_penetrance
from .structure import (
    compute_buried_surface_area,
    compute_interface_packing,
    compute_sasa,
    compute_shape_complementarity_approx,
    count_interface_contacts,
    count_interface_hbonds,
    count_interface_salt_bridges,
    find_hotspot_residues,
    find_interface_residues,
    get_chain_sequence,
    parse_structure,
)

logger = logging.getLogger(__name__)


def estimate_desolvation_energy(bsa: float, hydrophobic_bsa: float) -> float:
    """Eisenberg–McLachlan desolvation estimate (kcal/mol)."""
    polar_bsa = bsa - hydrophobic_bsa
    return round(-0.0278 * hydrophobic_bsa + 0.0160 * polar_bsa, 2)


def analyze_complex(
    pdb_path: str,
    binder_chain: str = "B",
    receptor_chain: str = "A",
    parental_pdb: Optional[str] = None,
    endogenous_pdb: Optional[str] = None,
    endogenous_binder_chain: Optional[str] = None,
    endogenous_receptor_chain: Optional[str] = None,
    foldx_bin: Optional[str] = None,
    foldx_repair: bool = False,
) -> Dict:
    """Analyse a single PDB complex and return a flat metrics dict."""
    result: Dict = {
        "filename": os.path.basename(pdb_path),
        "binder_chain": binder_chain,
        "receptor_chain": receptor_chain,
    }

    structure, model = parse_structure(pdb_path, "cx")

    chain_ids = [c.id for c in model]
    if binder_chain not in chain_ids or receptor_chain not in chain_ids:
        logger.error(
            "Chains %s/%s not in %s (available: %s)",
            binder_chain, receptor_chain, pdb_path, chain_ids,
        )
        return result

    # ---- Binder sequence ----
    binder_seq = get_chain_sequence(model, binder_chain)

    # ---- Rosetta metrics ----
    ros = compute_rosetta_metrics(pdb_path, binder_chain, receptor_chain)
    result["I_sc"] = ros.get("I_sc", float("nan"))
    result["dG_separated"] = ros.get("dG_separated", float("nan"))
    result["dSASA_int"] = ros.get("dSASA_int", float("nan"))
    sc_rosetta = ros.get("shape_complementarity", float("nan"))

    # dI_sc vs endogenous
    result["dI_sc"] = float("nan")
    if endogenous_pdb and os.path.exists(endogenous_pdb) and HAVE_PYROSETTA:
        ec_b = endogenous_binder_chain or binder_chain
        ec_r = endogenous_receptor_chain or receptor_chain
        endo = compute_rosetta_metrics(endogenous_pdb, ec_b, ec_r)
        if "I_sc" in endo:
            result["dI_sc"] = result["I_sc"] - endo["I_sc"]

    # ddG vs parental
    result["ddG_vs_parental"] = float("nan")
    if parental_pdb and os.path.exists(parental_pdb):
        result["ddG_vs_parental"] = compute_ddG_vs_parental(
            pdb_path, parental_pdb, f"{receptor_chain}_{binder_chain}"
        )

    # ---- FoldX metrics (fills gaps when PyRosetta unavailable) ----
    foldx_path = find_foldx_binary(foldx_bin)
    if foldx_path:
        logger.info("  FoldX …")
        fxm = analyse_complex_foldx(
            pdb_path, receptor_chain, binder_chain, foldx_path, repair=foldx_repair
        )
        if fxm:
            # Fill Rosetta-equivalent columns if still NaN
            if math.isnan(result.get("I_sc", float("nan"))):
                result["I_sc"] = fxm.get("foldx_interaction_energy", float("nan"))
            if math.isnan(result.get("dG_separated", float("nan"))):
                result["dG_separated"] = fxm.get("foldx_interaction_energy", float("nan"))
            if math.isnan(result.get("dSASA_int", float("nan"))):
                result["dSASA_int"] = fxm.get("foldx_desolvation", float("nan"))

            # dI_sc via FoldX if still NaN
            if math.isnan(result.get("dI_sc", float("nan"))) and \
               endogenous_pdb and os.path.exists(endogenous_pdb):
                ec_b = endogenous_binder_chain or binder_chain
                ec_r = endogenous_receptor_chain or receptor_chain
                endo_fx = analyse_complex_foldx(
                    endogenous_pdb, ec_r, ec_b, foldx_path, repair=foldx_repair
                )
                endo_ie = endo_fx.get("foldx_interaction_energy", float("nan"))
                my_ie = fxm.get("foldx_interaction_energy", float("nan"))
                try:
                    result["dI_sc"] = round(my_ie - endo_ie, 4)
                except (TypeError, ValueError):
                    pass

            # ddG via FoldX if still NaN
            if math.isnan(result.get("ddG_vs_parental", float("nan"))) and \
               parental_pdb and os.path.exists(parental_pdb):
                result["ddG_vs_parental"] = compute_foldx_ddG(
                    pdb_path, parental_pdb, receptor_chain, binder_chain, foldx_path
                )

            # Store all FoldX decomposed terms
            for k, v in fxm.items():
                result[k] = v

    # ---- PRODIGY binding affinity (open-source, no license) ----
    if HAVE_PRODIGY:
        logger.info("  PRODIGY …")
        pdg = compute_prodigy_metrics(pdb_path, receptor_chain, binder_chain)
        if pdg:
            # Fill dG if still NaN (Rosetta and FoldX both absent)
            if math.isnan(result.get("dG_separated", float("nan"))):
                result["dG_separated"] = pdg.get("prodigy_dG", float("nan"))
            if math.isnan(result.get("I_sc", float("nan"))):
                result["I_sc"] = pdg.get("prodigy_dG", float("nan"))

            # dI_sc via PRODIGY
            if math.isnan(result.get("dI_sc", float("nan"))) and \
               endogenous_pdb and os.path.exists(endogenous_pdb):
                ec_b = endogenous_binder_chain or binder_chain
                ec_r = endogenous_receptor_chain or receptor_chain
                endo_pdg = compute_prodigy_metrics(endogenous_pdb, ec_r, ec_b)
                endo_dg = endo_pdg.get("prodigy_dG", float("nan"))
                my_dg = pdg.get("prodigy_dG", float("nan"))
                try:
                    result["dI_sc"] = round(my_dg - endo_dg, 4)
                except (TypeError, ValueError):
                    pass

            # ddG via PRODIGY
            if math.isnan(result.get("ddG_vs_parental", float("nan"))) and \
               parental_pdb and os.path.exists(parental_pdb):
                result["ddG_vs_parental"] = compute_prodigy_ddG(
                    pdb_path, parental_pdb, receptor_chain, binder_chain
                )

            # Store all PRODIGY metrics
            for k, v in pdg.items():
                result[k] = v

    # ---- BSA ----
    logger.info("  BSA …")
    bsa, hydro_bsa = compute_buried_surface_area(pdb_path, receptor_chain, binder_chain)
    result["buried_surface_area"] = round(bsa, 1)
    result["hydrophobic_interface_area"] = round(hydro_bsa, 1)
    result["desolvation_energy"] = estimate_desolvation_energy(bsa, hydro_bsa)

    # ---- Shape complementarity ----
    if math.isnan(sc_rosetta):
        sc_rosetta = compute_shape_complementarity_approx(
            model, receptor_chain, binder_chain
        )
    result["shape_complementarity"] = (
        round(sc_rosetta, 3) if not math.isnan(sc_rosetta) else float("nan")
    )

    # ---- Contacts, packing, H-bonds, salt bridges ----
    result["n_interface_contacts"] = count_interface_contacts(
        model, receptor_chain, binder_chain
    )
    result["interface_packing_density"] = round(
        compute_interface_packing(model, receptor_chain, binder_chain), 2
    )

    if ros.get("n_hbonds_rosetta") is not None:
        result["n_hbonds"] = ros["n_hbonds_rosetta"]
    else:
        result["n_hbonds"] = count_interface_hbonds(model, receptor_chain, binder_chain)

    result["n_salt_bridges"] = count_interface_salt_bridges(
        model, receptor_chain, binder_chain
    )

    # Interface residue counts
    iface_b, iface_r = find_interface_residues(model, binder_chain, receptor_chain)
    result["n_interface_residues_binder"] = len(iface_b)
    result["n_interface_residues_receptor"] = len(iface_r)

    # ---- Hotspot residues ----
    result["hotspot_residues_binder"] = find_hotspot_residues(
        model, binder_chain, receptor_chain
    )

    # ---- SASA ----
    sasa = compute_sasa(pdb_path)
    result["total_sasa"] = round(sasa.get("total", 0), 1)
    result["binder_sasa"] = round(sasa.get(binder_chain, 0), 1)
    result["receptor_sasa"] = round(sasa.get(receptor_chain, 0), 1)

    # ---- Composition ----
    result.update(compute_composition(binder_seq))

    # ---- Secondary structure ----
    logger.info("  SS …")
    ss = compute_secondary_structure(pdb_path, model, binder_chain)
    ss_total = max(ss["total"], 1)
    result["n_helix_residues"] = ss["helix"]
    result["pct_helix"] = round(100 * ss["helix"] / ss_total, 1)
    result["n_sheet_residues"] = ss["sheet"]
    result["pct_sheet"] = round(100 * ss["sheet"] / ss_total, 1)
    result["n_loop_residues"] = ss["loop"]
    result["pct_loop"] = round(100 * ss["loop"] / ss_total, 1)

    # ---- Folds ----
    result["detected_folds"] = detect_folds(ss, binder_seq)

    # ---- Amphipathicity ----
    mean_uH, max_uH = compute_hydrophobic_moment(binder_seq)
    result["mean_hydrophobic_moment"] = round(mean_uH, 4)
    result["max_hydrophobic_moment"] = round(max_uH, 4)

    # ---- Cell penetrance ----
    result["predicted_cell_penetrance_score"] = predict_cell_penetrance(
        binder_seq,
        result["net_charge"],
        mean_uH,
        result["pct_hydrophobic"] / 100,
    )

    return result
