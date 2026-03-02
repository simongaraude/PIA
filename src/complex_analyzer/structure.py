"""BioPython structure helper functions."""

from __future__ import annotations

import os
import tempfile
from typing import Set, Tuple

import numpy as np
from Bio.PDB import PDBIO, NeighborSearch, PDBParser, Select, ShrakeRupley
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

try:
    from Bio.PDB.Polypeptide import three_to_one
except ImportError:
    from Bio.PDB.Polypeptide import protein_letters_3to1
    def three_to_one(resname: str) -> str:
        return protein_letters_3to1[resname]

from .constants import (
    ACCEPTOR_ATOMS,
    CONTACT_DISTANCE,
    DONOR_ATOMS,
    HBOND_DIST,
    HYDROPHOBIC_AA,
    INTERFACE_CB_DIST,
    NEG_SALT_ATOMS,
    POS_SALT_ATOMS,
    SALT_BRIDGE_DIST,
)

# ---------------------------------------------------------------------------
# Smart parser for PDB / CIF
# ---------------------------------------------------------------------------

def parse_structure(pdb_path: str, name: str = "s"):
    """Parse PDB or CIF file and return (structure, first_model)."""
    ext = os.path.splitext(pdb_path)[1].lower()
    if ext in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    struct = parser.get_structure(name, pdb_path)
    model = next(struct.get_models())
    return struct, model


# ---------------------------------------------------------------------------
# Basic accessors
# ---------------------------------------------------------------------------

def get_chain_sequence(model, chain_id: str) -> str:
    """Return 1-letter amino-acid sequence for *chain_id*."""
    seq: list[str] = []
    for residue in model[chain_id].get_residues():
        if is_aa(residue, standard=True):
            try:
                seq.append(three_to_one(residue.get_resname()))
            except KeyError:
                seq.append("X")
    return "".join(seq)


def get_chain_residues(model, chain_id: str):
    """Return list of standard amino-acid residues for *chain_id*."""
    return [r for r in model[chain_id].get_residues() if is_aa(r, standard=True)]


def _rep_atom(res):
    """Return CB (or CA for Gly) as representative atom."""
    if "CB" in res:
        return res["CB"]
    if "CA" in res:
        return res["CA"]
    return None


# ---------------------------------------------------------------------------
# Interface detection
# ---------------------------------------------------------------------------

def find_interface_residues(
    model, chain_a: str, chain_b: str, cutoff: float = INTERFACE_CB_DIST
) -> Tuple[Set, Set]:
    """Find interface residues using a CB–CB distance cutoff."""
    res_a = get_chain_residues(model, chain_a)
    res_b = get_chain_residues(model, chain_b)
    iface_a: Set = set()
    iface_b: Set = set()
    for ra in res_a:
        aa = _rep_atom(ra)
        if aa is None:
            continue
        for rb in res_b:
            ab = _rep_atom(rb)
            if ab is None:
                continue
            if aa - ab < cutoff:
                iface_a.add(ra)
                iface_b.add(rb)
    return iface_a, iface_b


# ---------------------------------------------------------------------------
# Contact counting
# ---------------------------------------------------------------------------

def count_interface_contacts(
    model, chain_a: str, chain_b: str, cutoff: float = CONTACT_DISTANCE
) -> int:
    """Count atom–atom contacts across the interface."""
    atoms_b = list(model[chain_b].get_atoms())
    ns = NeighborSearch(atoms_b)
    contacts = 0
    for atom in model[chain_a].get_atoms():
        contacts += len(ns.search(atom.get_vector().get_array(), cutoff))
    return contacts


# ---------------------------------------------------------------------------
# Hydrogen bonds (geometric criterion)
# ---------------------------------------------------------------------------

def count_interface_hbonds(
    model, chain_a: str, chain_b: str, dist_cutoff: float = HBOND_DIST
) -> int:
    """Estimate H-bonds across the interface (geometric donor–acceptor)."""
    atoms_b = list(model[chain_b].get_atoms())
    ns_b = NeighborSearch(atoms_b)
    hbonds = 0
    for atom_a in model[chain_a].get_atoms():
        name_a = atom_a.get_name().strip()
        is_donor = name_a in DONOR_ATOMS
        is_acceptor = name_a in ACCEPTOR_ATOMS
        if not (is_donor or is_acceptor):
            continue
        for atom_b in ns_b.search(atom_a.get_vector().get_array(), dist_cutoff):
            name_b = atom_b.get_name().strip()
            if is_donor and name_b in ACCEPTOR_ATOMS:
                hbonds += 1
            elif is_acceptor and name_b in DONOR_ATOMS:
                hbonds += 1
    return hbonds


# ---------------------------------------------------------------------------
# Salt bridges
# ---------------------------------------------------------------------------

def count_interface_salt_bridges(
    model, chain_a: str, chain_b: str, cutoff: float = SALT_BRIDGE_DIST
) -> int:
    """Count salt bridges across the interface."""
    def _charged_atoms(chain_id):
        pos, neg = [], []
        for res in model[chain_id].get_residues():
            if not is_aa(res, standard=True):
                continue
            rn = res.get_resname()
            for atom in res:
                n = atom.get_name().strip()
                if rn in ("LYS", "ARG", "HIS") and n in POS_SALT_ATOMS:
                    pos.append(atom)
                elif rn in ("ASP", "GLU") and n in NEG_SALT_ATOMS:
                    neg.append(atom)
        return pos, neg

    pos_a, neg_a = _charged_atoms(chain_a)
    pos_b, neg_b = _charged_atoms(chain_b)
    sb = 0
    for pa in pos_a:
        for nb in neg_b:
            if pa - nb < cutoff:
                sb += 1
    for pb in pos_b:
        for na in neg_a:
            if pb - na < cutoff:
                sb += 1
    return sb


# ---------------------------------------------------------------------------
# SASA helpers
# ---------------------------------------------------------------------------

class _ChainSelect(Select):
    """PDBIO selector that keeps only one chain."""
    def __init__(self, chain_id: str):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id


def compute_sasa(pdb_path: str) -> dict:
    """Per-chain and total SASA via BioPython Shrake–Rupley."""
    struct, model = parse_structure(pdb_path, "s")
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    result: dict = {}
    total = 0.0
    for chain in model:
        chain_sasa = sum(r.sasa for r in chain if is_aa(r, standard=True))
        result[chain.id] = chain_sasa
        total += chain_sasa
    result["total"] = total
    return result


def compute_buried_surface_area(
    pdb_path: str, chain_a: str, chain_b: str
) -> Tuple[float, float]:
    """
    BSA = SASA_A(alone) + SASA_B(alone) − SASA_AB(complex).

    Returns ``(bsa, hydrophobic_bsa)`` in Å².
    """
    struct, model = parse_structure(pdb_path, "cx")
    sr = ShrakeRupley()

    # Complex SASA
    sr.compute(model, level="A")
    sasa_complex = sum(a.sasa for a in model.get_atoms())

    io = PDBIO()

    def _isolated_sasa(chain_id: str) -> float:
        io.set_structure(struct)
        # Always write isolated chains as PDB (PDBIO outputs PDB format)
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        try:
            io.save(tmp.name, _ChainSelect(chain_id))
            iso_parser = PDBParser(QUIET=True)
            s = iso_parser.get_structure("iso", tmp.name)
            iso_model = next(s.get_models())
            sr.compute(iso_model, level="A")
            return sum(a.sasa for a in iso_model.get_atoms())
        finally:
            os.unlink(tmp.name)

    sasa_a_alone = _isolated_sasa(chain_a)
    sasa_b_alone = _isolated_sasa(chain_b)
    bsa = sasa_a_alone + sasa_b_alone - sasa_complex

    # Hydrophobic fraction at interface → hydrophobic BSA
    iface_a, iface_b = find_interface_residues(model, chain_a, chain_b, cutoff=8.0)
    n_iface = len(iface_a) + len(iface_b)
    n_hydrophobic = 0
    for r in list(iface_a) + list(iface_b):
        try:
            if three_to_one(r.get_resname()) in HYDROPHOBIC_AA:
                n_hydrophobic += 1
        except KeyError:
            pass
    frac = n_hydrophobic / max(n_iface, 1)
    return bsa, bsa * frac


# ---------------------------------------------------------------------------
# Shape complementarity (approximate)
# ---------------------------------------------------------------------------

def compute_shape_complementarity_approx(
    model, chain_a: str, chain_b: str
) -> float:
    """
    Approximate Sc via local surface-normal dot products.

    For publication-quality Sc, use Rosetta InterfaceAnalyzerMover.
    """
    atoms_a = [a for r in get_chain_residues(model, chain_a) for a in r.get_atoms()]
    atoms_b = [a for r in get_chain_residues(model, chain_b) for a in r.get_atoms()]
    ns_a = NeighborSearch(atoms_a)
    ns_b = NeighborSearch(atoms_b)

    pairs = []
    for a in atoms_a:
        for b in ns_b.search(a.get_vector().get_array(), 5.0):
            pairs.append((a, b))
    if len(pairs) < 10:
        return float("nan")

    def _normal(atom, ns, k=6):
        coord = atom.get_vector().get_array()
        nbrs = ns.search(coord, 6.0)
        if len(nbrs) < 3:
            return None
        centroid = np.mean([n.get_vector().get_array() for n in nbrs[:k]], axis=0)
        v = coord - centroid
        norm = np.linalg.norm(v)
        return v / norm if norm > 1e-6 else None

    dots = []
    for a, b in pairs[:500]:
        na = _normal(a, ns_a)
        nb = _normal(b, ns_b)
        if na is not None and nb is not None:
            dots.append(max(0.0, -np.dot(na, nb)))

    return float(np.median(dots)) if dots else float("nan")


# ---------------------------------------------------------------------------
# Interface packing density
# ---------------------------------------------------------------------------

def compute_interface_packing(model, chain_a: str, chain_b: str) -> float:
    """Mean atom neighbors within 4 Å for interface atoms."""
    iface_a, iface_b = find_interface_residues(model, chain_a, chain_b)
    iface_atoms = [a for r in list(iface_a) + list(iface_b) for a in r.get_atoms()]
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)
    densities = [len(ns.search(a.get_vector().get_array(), 4.0)) for a in iface_atoms]
    return float(np.mean(densities)) if densities else float("nan")


# ---------------------------------------------------------------------------
# Hotspot residues (contact-count proxy)
# ---------------------------------------------------------------------------

def find_hotspot_residues(
    model, binder_chain: str, receptor_chain: str, n_top: int = 10
) -> str:
    """Top-N binder residues by cross-interface atom contacts."""
    iface_binder, _ = find_interface_residues(model, binder_chain, receptor_chain)
    atoms_receptor = list(model[receptor_chain].get_atoms())
    ns = NeighborSearch(atoms_receptor)
    scores = []
    for res in iface_binder:
        contacts = sum(
            len(ns.search(a.get_vector().get_array(), CONTACT_DISTANCE))
            for a in res
        )
        try:
            aa = three_to_one(res.get_resname())
        except KeyError:
            aa = "X"
        scores.append((f"{aa}{res.get_id()[1]}", contacts))
    scores.sort(key=lambda x: x[1], reverse=True)
    return ", ".join(f"{r}({c}ct)" for r, c in scores[:n_top])
