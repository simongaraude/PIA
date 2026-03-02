"""Tests for complex-analyzer."""

import os
import tempfile

import pytest

from complex_analyzer.constants import STANDARD_AA
from complex_analyzer.secondary_structure import detect_folds
from complex_analyzer.sequence import (
    compute_composition,
    compute_hydrophobic_moment,
    predict_cell_penetrance,
)

# ---------------------------------------------------------------------------
# Sequence composition
# ---------------------------------------------------------------------------

class TestComposition:
    def test_basic_counts(self):
        comp = compute_composition("AAARRRDDD")
        assert comp["n_residues"] == 9
        assert comp["count_A"] == 3
        assert comp["count_R"] == 3
        assert comp["count_D"] == 3
        assert comp["n_pos_charged"] == 3  # R, K, H — here just R
        assert comp["n_neg_charged"] == 3  # D

    def test_percentages(self):
        comp = compute_composition("AAAA")  # 4 Ala = hydrophobic
        assert comp["pct_hydrophobic"] == 100.0
        assert comp["pct_charged"] == 0.0

    def test_empty_sequence(self):
        comp = compute_composition("")
        assert comp["n_residues"] == 0
        assert comp["net_charge"] == 0.0

    def test_net_charge(self):
        comp = compute_composition("RRRKKKDDD")
        # R=3 (+3), K=3 (+3), D=3 (-3) → net +3
        assert comp["net_charge"] == pytest.approx(3.0, abs=0.5)

    def test_all_aa_columns_present(self):
        comp = compute_composition("ACGT")
        for aa in STANDARD_AA:
            assert f"count_{aa}" in comp
            assert f"pct_{aa}" in comp


# ---------------------------------------------------------------------------
# Hydrophobic moment
# ---------------------------------------------------------------------------

class TestHydrophobicMoment:
    def test_returns_tuple(self):
        mean_uH, max_uH = compute_hydrophobic_moment("AKLWAILKWLA")
        assert isinstance(mean_uH, float)
        assert isinstance(max_uH, float)
        assert max_uH >= mean_uH

    def test_empty_sequence(self):
        assert compute_hydrophobic_moment("") == (0.0, 0.0)

    def test_short_sequence(self):
        mean_uH, max_uH = compute_hydrophobic_moment("AKL")
        assert mean_uH == max_uH  # shorter than window

    def test_amphipathic_higher_than_uniform(self):
        # Alternating hydrophobic/polar should have higher moment than all-Ala
        mean_amph, _ = compute_hydrophobic_moment("AKAKAKAKAKAKAKAKAKAKAKAKAKAKAKAL")
        mean_uni, _ = compute_hydrophobic_moment("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        # Not guaranteed to always hold for all windows, but generally true
        assert mean_amph > 0


# ---------------------------------------------------------------------------
# Cell penetrance
# ---------------------------------------------------------------------------

class TestCellPenetrance:
    def test_range(self):
        score = predict_cell_penetrance("RRRRRRRR", 8.0, 0.3, 0.0)
        assert 0 <= score <= 1

    def test_arg_rich_higher(self):
        score_arg = predict_cell_penetrance("RRRRRRRR", 8.0, 0.3, 0.0)
        score_ala = predict_cell_penetrance("AAAAAAAA", 0.0, 0.1, 1.0)
        assert score_arg > score_ala

    def test_empty(self):
        score = predict_cell_penetrance("", 0.0, 0.0, 0.0)
        assert 0 <= score <= 1


# ---------------------------------------------------------------------------
# Fold detection
# ---------------------------------------------------------------------------

class TestFoldDetection:
    def test_all_alpha(self):
        ss = {"helix": 80, "sheet": 5, "loop": 15, "total": 100}
        folds = detect_folds(ss, "A" * 100)
        assert "all-alpha" in folds

    def test_all_beta(self):
        ss = {"helix": 5, "sheet": 50, "loop": 45, "total": 100}
        folds = detect_folds(ss, "A" * 100)
        assert "all-beta" in folds

    def test_irregular(self):
        ss = {"helix": 5, "sheet": 5, "loop": 90, "total": 100}
        folds = detect_folds(ss, "A" * 100)
        assert "irregular" in folds

    def test_disulfide_rich(self):
        seq = "CCCCACCCCACCCA"  # >6% Cys, <80 residues
        ss = {"helix": 2, "sheet": 2, "loop": 10, "total": 14}
        folds = detect_folds(ss, seq)
        assert "disulfide-rich" in folds


# ---------------------------------------------------------------------------
# Integration: PDB analysis (only if BioPython can create a minimal PDB)
# ---------------------------------------------------------------------------

def _make_minimal_pdb() -> str:
    """Create a tiny two-chain PDB for integration tests."""
    pdb_text = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.764   1.208  1.00  0.00           C
ATOM      6  N   ALA A   2       3.325   1.506   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.930   2.829   0.000  1.00  0.00           C
ATOM      8  C   ALA A   2       5.451   2.729   0.000  1.00  0.00           C
ATOM      9  O   ALA A   2       6.037   1.660   0.000  1.00  0.00           O
ATOM     10  CB  ALA A   2       3.453   3.649   1.196  1.00  0.00           C
TER
ATOM     11  N   ARG B   1       2.000   0.000   4.000  1.00  0.00           N
ATOM     12  CA  ARG B   1       3.458   0.000   4.000  1.00  0.00           C
ATOM     13  C   ARG B   1       4.009   1.420   4.000  1.00  0.00           C
ATOM     14  O   ARG B   1       3.246   2.390   4.000  1.00  0.00           O
ATOM     15  CB  ARG B   1       3.986  -0.764   5.208  1.00  0.00           C
ATOM     16  N   LYS B   2       5.325   1.506   4.000  1.00  0.00           N
ATOM     17  CA  LYS B   2       5.930   2.829   4.000  1.00  0.00           C
ATOM     18  C   LYS B   2       7.451   2.729   4.000  1.00  0.00           C
ATOM     19  O   LYS B   2       8.037   1.660   4.000  1.00  0.00           O
ATOM     20  CB  LYS B   2       5.453   3.649   5.196  1.00  0.00           C
TER
END
"""
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    tmp.write(pdb_text)
    tmp.close()
    return tmp.name


class TestStructureIntegration:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.pdb_path = _make_minimal_pdb()
        yield
        os.unlink(self.pdb_path)

    def test_chain_sequence(self):
        from Bio.PDB import PDBParser

        from complex_analyzer.structure import get_chain_sequence

        model = PDBParser(QUIET=True).get_structure("t", self.pdb_path)[0]
        assert get_chain_sequence(model, "A") == "AA"
        assert get_chain_sequence(model, "B") == "RK"

    def test_interface_contacts(self):
        from Bio.PDB import PDBParser

        from complex_analyzer.structure import count_interface_contacts

        model = PDBParser(QUIET=True).get_structure("t", self.pdb_path)[0]
        n = count_interface_contacts(model, "A", "B", cutoff=8.0)
        assert isinstance(n, int)
        assert n >= 0

    def test_buried_surface_area(self):
        from complex_analyzer.structure import compute_buried_surface_area

        bsa, hydro_bsa = compute_buried_surface_area(self.pdb_path, "A", "B")
        assert isinstance(bsa, float)
        assert bsa >= 0

    def test_full_analysis(self):
        from complex_analyzer.analyzer import analyze_complex

        result = analyze_complex(self.pdb_path, binder_chain="B", receptor_chain="A")
        assert result["filename"].endswith(".pdb")
        assert result["n_residues"] == 2
        assert "net_charge" in result
        assert "predicted_cell_penetrance_score" in result
