#!/usr/bin/env python3

import sys
import json
from collections import defaultdict
from Bio.PDB import MMCIFParser


if len(sys.argv) != 3:
    sys.exit(
        f"Usage: {sys.argv[0]} <confidences.json> <model.cif>"
    )

json_file = sys.argv[1]
cif_file = sys.argv[2]

# ------------------------------------------------------------
# Read AF3 confidence JSON
# ------------------------------------------------------------

with open(json_file) as f:
    data = json.load(f)

atom_chain_ids = data["atom_chain_ids"]
atom_plddts = data["atom_plddts"]

# ------------------------------------------------------------
# Read CIF
# ------------------------------------------------------------

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("model", cif_file)

# ------------------------------------------------------------
# Get atoms from CIF
# ------------------------------------------------------------

atoms = list(structure.get_atoms())

if len(atoms) != len(atom_plddts):
    sys.exit(
        f"ERROR: CIF contains {len(atoms)} atoms but JSON "
        f"contains {len(atom_plddts)} atom pLDDTs"
    )

if len(atoms) != len(atom_chain_ids):
    sys.exit(
        f"ERROR: CIF contains {len(atoms)} atoms but JSON "
        f"contains {len(atom_chain_ids)} atom chain IDs"
    )

# ------------------------------------------------------------
# Collect pLDDT values for each residue in chain A
# ------------------------------------------------------------

residue_plddts = defaultdict(list)
residue_aa = {}

for atom, chain_id, plddt in zip(
    atoms,
    atom_chain_ids,
    atom_plddts
):

    # Only analyse probe chain A
    if chain_id != "A":
        continue

    residue = atom.get_parent()

    hetflag, resseq, icode = residue.id

    # Ignore waters / hetero residues
    if hetflag != " ":
        continue

    residue_plddts[resseq].append(plddt)

    # Keep native 3-letter amino-acid code from CIF
    residue_aa[resseq] = residue.get_resname()

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

print("residue\taa\tpLDDT")

for residue in sorted(residue_plddts):

    values = residue_plddts[residue]
    mean_plddt = sum(values) / len(values)

    print(
        f"{residue}\t"
        f"{residue_aa[residue]}\t"
        f"{mean_plddt:.2f}"
    )