# PIA (Protein Interface Analyzer)

Batch structural analysis of protein-protein complexes from computational design pipelines. Reads PDB/CIF files and produces a single CSV with 80+ biophysical, energetic, and compositional metrics per complex.

Built for high-throughput scoring of binder design campaigns (BindCraft, BoltzGen, RFdiffusion, etc.) on AWS ParallelCluster with Slurm. Supports automatic Google Drive sync via rclone.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Usage](#usage)
4. [AWS ParallelCluster Deployment](#aws-parallelcluster-deployment)
5. [Metrics Reference](#metrics-reference)
   - [Binding Affinity (PRODIGY)](#binding-affinity-prodigy)
   - [Interface Geometry](#interface-geometry)
   - [Energetics](#energetics)
   - [Binder Characterisation](#binder-characterisation)
   - [Sequence Composition](#sequence-composition)
   - [Secondary Structure and Fold Classification](#secondary-structure-and-fold-classification)
   - [Amphipathicity](#amphipathicity)
   - [Cell Penetrance](#cell-penetrance)
6. [Dependency Tiers](#dependency-tiers)
7. [Project Structure](#project-structure)
8. [Caveats and Limitations](#caveats-and-limitations)
9. [References](#references)
10. [License](#license)

---

## Quick Start

```bash
pip install biopython pandas numpy freesasa prodigy-prot
pip install -e .
complex-analyzer --input_dir ./my_structures/ -o results.csv
```

---

## Installation

```bash
git clone https://github.com/YOUR_USER/complex-analyzer.git
cd complex-analyzer

pip install biopython pandas numpy freesasa prodigy-prot
pip install -e .

# Optional: DSSP for dictionary-based secondary structure assignment
sudo apt install dssp
```

### Development

```bash
pip install -e ".[dev]"
pytest -v
ruff check src/ tests/
```

---

## Usage

### Command line

```bash
# Analyse all PDB/CIF files in a directory
complex-analyzer --input_dir ./structures/ -o metrics.csv

# Analyse specific files
complex-analyzer --input_files complex1.pdb complex2.cif -o metrics.csv

# Custom chain IDs (default: A = receptor, B = binder)
complex-analyzer --input_dir ./structures/ -r A -b B

# With reference structures for differential scoring
complex-analyzer --input_dir ./structures/ \
    --endogenous_pdb reference/endogenous.pdb \
    --endogenous_binder_chain C
```

### As a Python library

```python
from complex_analyzer.analyzer import analyze_complex

metrics = analyze_complex(
    "complex.pdb",
    binder_chain="B",
    receptor_chain="A",
)

print(metrics["prodigy_dG"])
print(metrics["buried_surface_area"])
print(metrics["predicted_cell_penetrance_score"])
```

---

## AWS ParallelCluster Deployment

These instructions assume an AWS ParallelCluster with Slurm, shared storage at `/shared/`, and a conda environment on shared storage.

### 1. Upload and install

From your local machine:

```bash
scp -i ~/.ssh/your_key.pem complex-analyzer.zip ubuntu@YOUR_HEADNODE_IP:/shared/
```

On the headnode:

```bash
cd /shared && unzip complex-analyzer.zip -d complex-analyzer
source /shared/miniconda3/etc/profile.d/conda.sh && conda activate your_env
bash /shared/complex-analyzer/scripts/install_analyzer.sh
```

### 2. Prepare structures

Place all PDB or CIF files into a single flat directory:

```bash
mkdir -p /shared/all_structures
find /path/to/outputs -name "*.cif" -exec cp {} /shared/all_structures/ \;
```

### 3. Submit the analysis job

```bash
/shared/complex-analyzer/scripts/submit_analysis.sh /shared/all_structures
```

### 4. Monitor

```bash
squeue
tail -f /shared/all_structures/analysis_job/slurm_<ID>.err
```

### 5. Results

Output CSV at `/shared/all_structures/complex_metrics.csv`, optionally synced to Google Drive via rclone.

### Configuring rclone for Google Drive

```bash
curl https://rclone.org/install.sh | sudo bash
rclone config   # follow interactive prompts for Google Drive
```

Update `GDRIVE_REMOTE` and `GDRIVE_FOLDER` in `scripts/submit_analysis.sh`.

---

## Metrics Reference

### Binding Affinity (PRODIGY)

Binding affinity prediction using PRODIGY (PROtein binDIng enerGY prediction; Xue et al., Bioinformatics 2016). PRODIGY is open-source (Apache 2.0), requires no license, and predicts binding free energy from the 3D structure of the complex using a linear model trained on intermolecular contact statistics and non-interacting surface (NIS) properties.

| Column | Description |
|--------|-------------|
| `prodigy_dG` | Predicted binding free energy (kcal/mol). More negative values indicate stronger predicted binding. |
| `prodigy_n_contacts` | Total number of intermolecular contacts at a 5.5 A distance cutoff between heavy atoms of interface residues. |
| `prodigy_cc` | Count of charged-charged intermolecular contacts (Arg, Lys, His, Asp, Glu on both sides). |
| `prodigy_cp` | Count of charged-polar contacts. |
| `prodigy_ca` | Count of charged-apolar contacts. |
| `prodigy_pp` | Count of polar-polar contacts. |
| `prodigy_ap` | Count of apolar-polar contacts. |
| `prodigy_aa` | Count of apolar-apolar contacts. |
| `prodigy_nis_apolar` | Percentage of apolar residues on the non-interacting surface (NIS). NIS residues are surface-exposed but not at the interface. This captures the hydrophobic context surrounding the binding site. |
| `prodigy_nis_charged` | Percentage of charged residues on the non-interacting surface. |

The PRODIGY linear model: dG = 0.09459 * n_cc + 0.01127 * n_aa + ... (contact type weights) + NIS terms. Benchmark performance: Pearson r = 0.73, RMSE = 1.89 kcal/mol against experimental binding affinities.

When PRODIGY is installed, `I_sc` and `dG_separated` are filled with `prodigy_dG` as a fallback if PyRosetta and FoldX are not available.

### Interface Geometry

| Column | Description |
|--------|-------------|
| `buried_surface_area` | Total buried surface area (A^2) upon complex formation, computed as (SASA_receptor + SASA_binder - SASA_complex) / 2. Uses FreeSASA if available, otherwise BioPython's Shrake-Rupley algorithm with probe radius 1.4 A. |
| `hydrophobic_interface_area` | Estimated hydrophobic contribution to BSA (A^2), computed by weighting buried residue areas by Kyte-Doolittle hydrophobicity. |
| `shape_complementarity` | Geometric shape complementarity of the interface (0 to 1). Without PyRosetta, estimated by computing surface normal vectors at interface atoms and measuring the alignment between opposing surface patches. A value of 1 indicates perfect geometric complementarity. With PyRosetta, uses Rosetta's implementation of the Lawrence and Colman (1993) Sc metric. |
| `n_interface_contacts` | Number of inter-chain atom pairs within 4.5 A. |
| `interface_packing_density` | Mean number of heavy-atom neighbours within 4.0 A for interface residues. Higher values indicate tighter packing. |
| `n_hbonds` | Number of inter-chain hydrogen bonds, detected geometrically: donor-acceptor distance < 3.5 A, D-H...A angle > 120 degrees. |
| `n_salt_bridges` | Number of inter-chain salt bridges: oppositely charged atom pairs (e.g., Lys NZ to Asp OD1) within 4.0 A. |
| `n_interface_residues_binder` | Number of binder residues with at least one atom within 4.5 A of the receptor. |
| `n_interface_residues_receptor` | Number of receptor residues with at least one atom within 4.5 A of the binder. |
| `hotspot_residues_binder` | Top 10 binder residues ranked by number of cross-interface atomic contacts. Format: residue name, chain, sequence position. |

### Energetics

| Column | Description |
|--------|-------------|
| `I_sc` | Interface score. From PyRosetta InterfaceAnalyzerMover if available; otherwise filled with `prodigy_dG`. |
| `dI_sc` | Differential interface score relative to an endogenous complex. Requires `--endogenous_pdb`. |
| `dG_separated` | Binding energy from separating the two chains. From PyRosetta (separate + repack protocol) or `prodigy_dG` as fallback. |
| `dSASA_int` | Change in solvent-accessible surface area upon binding. From PyRosetta. |
| `desolvation_energy` | Estimated desolvation penalty (kcal/mol) using the Eisenberg-McLachlan solvation model. Computed as the sum over buried atoms of (atomic solvation parameter * change in SASA). Positive values indicate an energetic cost of desolvation; lower values are more favorable. |

### Binder Characterisation

| Column | Description |
|--------|-------------|
| `n_residues` | Number of residues in the binder chain. |
| `total_sasa` | Total SASA of the complex (A^2). |
| `binder_sasa` | SASA of the binder chain in isolation (A^2). |
| `receptor_sasa` | SASA of the receptor chain in isolation (A^2). |
| `net_charge` | Net charge of the binder at pH 7. Computed as: +1 per Arg and Lys, -1 per Asp and Glu, +0.1 per His (partial protonation). |

### Sequence Composition

All composition metrics are computed on the binder chain only.

| Column | Description |
|--------|-------------|
| `n_charged`, `pct_charged` | Count and percentage of charged residues (R, K, H, D, E). |
| `n_pos_charged`, `pct_pos_charged` | Count and percentage of positively charged residues (R, K, H). |
| `n_neg_charged`, `pct_neg_charged` | Count and percentage of negatively charged residues (D, E). |
| `n_hydrophobic`, `pct_hydrophobic` | Count and percentage of hydrophobic residues (A, I, L, M, F, W, V, P). |
| `n_polar`, `pct_polar` | Count and percentage of polar residues (S, T, N, Q, Y, C). |
| `count_X`, `pct_X` | Count and percentage for each of the 20 standard amino acids (A through Y). |

### Secondary Structure and Fold Classification

| Column | Description |
|--------|-------------|
| `n_helix_residues`, `pct_helix` | Count and percentage of residues in alpha-helical conformation. |
| `n_sheet_residues`, `pct_sheet` | Count and percentage of residues in beta-sheet conformation. |
| `n_loop_residues`, `pct_loop` | Count and percentage of residues in loop/coil conformation. |
| `detected_folds` | Fold classification based on secondary structure content. |

Secondary structure is assigned by DSSP (Kabsch and Sander, 1983) when `mkdssp` is installed. DSSP eight-state codes are reduced to three states: H, G, I mapped to helix; E, B mapped to sheet; all others mapped to loop.

When DSSP is unavailable, a CA-distance heuristic is used: consecutive CA-CA distances of 5.2-5.6 A suggest beta-strand geometry; CA(i) to CA(i+4) distances of 5.0-6.5 A suggest helical geometry.

Fold classification rules applied to the binder chain:
- All-alpha: helix content > 40% and sheet content < 10%
- All-beta: sheet content > 30% and helix content < 10%
- Alpha+beta: both helix > 15% and sheet > 15%
- Coiled-coil: helix content > 70% and binder length > 40 residues
- Disulfide-rich: cysteine count >= 4 and binder length < 80 residues

### Amphipathicity

Amphipathicity is quantified by the hydrophobic moment (Eisenberg et al., 1982), which measures the periodicity of hydrophobic and hydrophilic residues along the sequence. A high hydrophobic moment indicates that one face of a helix (or strand) is hydrophobic while the opposite face is hydrophilic -- a hallmark of membrane-interacting or protein-binding surfaces.

| Column | Description |
|--------|-------------|
| `mean_hydrophobic_moment` | Mean hydrophobic moment across all sliding windows. |
| `max_hydrophobic_moment` | Maximum hydrophobic moment across all sliding windows. |

The calculation uses a sliding window of 11 residues (default) with a helical rotation angle of 100 degrees (appropriate for alpha-helices). For each window of residues r_1 through r_n, the hydrophobic moment mu_H is:

```
mu_H = (1/n) * sqrt( [sum_i H(r_i) * sin(i * delta)]^2 + [sum_i H(r_i) * cos(i * delta)]^2 )
```

where H(r_i) is the Eisenberg consensus hydrophobicity of residue r_i, delta = 100 degrees is the angular separation between consecutive residues in an alpha-helix, and n is the window length. The Eisenberg consensus scale ranges from -2.53 (Arg, most hydrophilic) to +1.38 (Ile, most hydrophobic).

The mean hydrophobic moment is the average of mu_H across all windows; the max is the single highest window value. Higher values (> 0.3-0.4) suggest amphipathic helical segments. Values near zero suggest either uniform hydrophobicity or uniform hydrophilicity.

### Cell Penetrance

| Column | Description |
|--------|-------------|
| `predicted_cell_penetrance_score` | Heuristic cell-penetrance score (0 to 1). |

This is a weighted heuristic inspired by features used in validated cell-penetrating peptide (CPP) predictors such as CPPpred (Holton et al., 2013), CellPPD (Gautam et al., 2013), and MLCPP (Manavalan et al., 2018). It is not a validated predictor and should be used only for rough ranking within a campaign, not for absolute predictions.

The score is computed as a weighted sum of six features, each normalised to [0, 1]:

1. **Positive charge (weight: 0.25)**: CPPs typically carry a net positive charge between +4 and +12. Score = 0.25 * min(net_charge / 8, 1). Capped at +8 net charge.

2. **Arginine content (weight: 0.20)**: Arginine-rich peptides are strongly associated with membrane translocation due to guanidinium-membrane interactions. Score = 0.20 * min(Arg_fraction / 0.25, 1). Saturates when Arg exceeds 25% of the sequence.

3. **Amphipathicity (weight: 0.20)**: Amphipathic helices can insert into membranes. Score = 0.20 * min(mean_hydrophobic_moment / 0.5, 1). Uses the mean hydrophobic moment computed above.

4. **Moderate hydrophobicity (weight: 0.15)**: CPPs are neither strongly hydrophobic nor strongly hydrophilic. The score peaks when mean Kyte-Doolittle hydrophobicity is near 0.5 and decreases with deviation. Score = 0.15 * max(0, 1 - |mean_KD - 0.5| / 3).

5. **Length penalty (weight: 0.10)**: Shorter peptides penetrate membranes more readily. Score = 0.10 for peptides of 30 residues or fewer, 0.05 for 31-60, 0.02 for 61-100, 0.005 for longer sequences.

6. **Aromatic anchoring (weight: 0.10)**: Trp and Phe residues anchor peptides at the membrane-water interface. Score = 0.10 * min(aromatic_fraction / 0.15, 1), where aromatic_fraction = (count_W + count_F) / n_residues.

The total score is clamped to [0, 1]. A score above 0.5 suggests the binder has multiple features associated with cell penetration, but experimental validation is required.

---

## Dependency Tiers

The tool degrades gracefully depending on what is installed.

| Tier | Packages | Metrics gained |
|------|----------|----------------|
| Minimum | biopython, pandas, numpy | All structural metrics with geometric approximations |
| Recommended | + freesasa, prodigy-prot | Accurate SASA; binding affinity (dG, Kd, contact decomposition) |
| Full (DSSP) | + dssp (mkdssp) | Dictionary-based secondary structure assignment |
| Full (PyRosetta) | + PyRosetta | I_sc, dG_separated, dSASA, ddG, Lawrence-Colman Sc |
| Full (FoldX) | + FoldX binary | Energy decomposition, interaction energy, van der Waals, electrostatics |

Priority cascade for binding energy columns: PyRosetta > FoldX > PRODIGY.

---

## Project Structure

```
complex-analyzer/
  src/complex_analyzer/
    __init__.py
    __main__.py
    cli.py                  # CLI entry point and column ordering
    analyzer.py             # Orchestrator
    structure.py            # Structural metrics (BSA, contacts, Sc, H-bonds)
    sequence.py             # Composition, amphipathicity, cell penetrance
    secondary_structure.py  # DSSP / CA-heuristic / fold classification
    rosetta_metrics.py      # Optional PyRosetta integration
    foldx_metrics.py        # Optional FoldX integration
    prodigy_metrics.py      # PRODIGY binding affinity
    constants.py            # Hydrophobicity scales, cutoffs, atom sets
  scripts/
    install_analyzer.sh     # One-time install on cluster
    submit_analysis.sh      # Slurm job submission
  tests/
    test_analyzer.py
  .github/workflows/
    ci.yml
  pyproject.toml
  requirements.txt
  LICENSE
```

---

## Caveats and Limitations

**Shape complementarity** -- Without PyRosetta, uses a surface-normal dot-product approximation. Suitable for ranking designs within a campaign but not for quantitative comparison across studies. For publication-quality Sc values, use Rosetta's InterfaceAnalyzerMover (Lawrence and Colman, J. Mol. Biol. 1993).

**PRODIGY binding affinity** -- A linear model trained on experimental binding affinities. Good for relative ranking within a set of designs targeting the same receptor. Not suitable for comparing absolute affinities across different targets. Performance degrades for complexes outside the training domain (e.g., very small peptides, intrinsically disordered regions).

**Cell penetrance score** -- A heuristic, not a validated predictor. For experimental planning, submit sequences to CPPpred (http://bioware.ucd.ie/~compass/biowareweb/Server_pages/cpppred.php), CellPPD (https://webs.iiitd.edu.in/raghava/cellppd/), or MLCPP (http://www.thegleelab.org/MLCPP/).

**Fold classification** -- Based solely on secondary structure content thresholds. Does not perform structural alignment. For proper fold assignment, use Foldseek (van Kempen et al., Nature Biotechnology 2024) or classify against CATH/SCOPe.

**Hotspot residues** -- Ranked by contact count, not by energetic contribution. For rigorous computational alanine scanning, use Rosetta ddg_monomer or FoldX PositionScan.

**Secondary structure without DSSP** -- The CA-distance heuristic is approximate. It tends to undercount short helices and may misassign residues at helix/sheet boundaries. Install `mkdssp` for accurate assignment.

**CIF/mmCIF support** -- All mmCIF files are read via BioPython's MMCIFParser. DSSP and PRODIGY CLI require PDB format, so CIF files are automatically converted on the fly using BioPython's PDBIO. This conversion may lose some metadata but preserves coordinates and chain assignments.

---

## License

MIT
