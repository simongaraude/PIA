"""Shared constants and lookup tables."""

# Distance cutoffs (Angstroms)
CONTACT_DISTANCE = 4.5
HBOND_DIST = 3.5
HBOND_ANGLE_MIN = 120.0
SALT_BRIDGE_DIST = 4.0
INTERFACE_CB_DIST = 10.0

# Amino acid classifications
HYDROPHOBIC_AA = set("AILMFWVP")
POLAR_AA = set("STNQYC")
POS_CHARGED_AA = set("RKH")
NEG_CHARGED_AA = set("DE")
CHARGED_AA = POS_CHARGED_AA | NEG_CHARGED_AA
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydrophobicity scale
KD_HYDROPHOBICITY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Eisenberg consensus hydrophobicity (for hydrophobic moment)
EISENBERG_HYDROPHOBICITY = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

# Net charge at pH 7 (approximate)
AA_CHARGE = {
    "R": 1, "K": 1, "H": 0.1,
    "D": -1, "E": -1,
}

# Atom name sets for H-bond / salt-bridge detection
DONOR_ATOMS = {
    "N", "NE", "NE1", "NE2", "ND1", "ND2",
    "NH1", "NH2", "NZ", "OG", "OG1", "OH", "SG",
}
ACCEPTOR_ATOMS = {
    "O", "OD1", "OD2", "OE1", "OE2",
    "OG", "OG1", "OH", "SD", "ND1", "NE2",
}
POS_SALT_ATOMS = {"NZ", "NH1", "NH2", "NE", "ND1", "NE2"}
NEG_SALT_ATOMS = {"OD1", "OD2", "OE1", "OE2"}
