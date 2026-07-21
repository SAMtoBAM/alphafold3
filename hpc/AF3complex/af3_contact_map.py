#!/usr/bin/env python3

import json
import argparse
import numpy as np
import csv

from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO


def main():

    parser = argparse.ArgumentParser(
        description="Map AlphaFold3 contact probabilities onto an mmCIF model."
    )

    parser.add_argument(
        "--json",
        required=True,
        help="AlphaFold3 full_data/confidences JSON"
    )

    parser.add_argument(
        "--cif",
        required=True,
        help="AlphaFold3 model mmCIF"
    )

    parser.add_argument(
        "--out_cif",
        default="interface_model.cif",
        help="Output mmCIF with interface scores stored as B-factors"
    )

    parser.add_argument(
        "--out_tsv",
        default="interface_scores.tsv",
        help="Output TSV of residue interface scores"
    )

    parser.add_argument(
        "--partner_threshold",
        type=float,
        default=0.5,
        help="Minimum contact probability for reporting interface partners"
    )

    args = parser.parse_args()

    ############################################################
    # Read JSON
    ############################################################

    print("Loading JSON...")

    with open(args.json) as f:
        data = json.load(f)

    contact_probs = np.array(data["contact_probs"])
    chains = np.array(data["token_chain_ids"])
    residues = np.array(data["token_res_ids"])

    print(f"Contact matrix: {contact_probs.shape}")

    ############################################################
    # Calculate residue interface scores
    ############################################################

    print("Calculating interface scores...")

    scores = {}

    unique_chains = np.unique(chains)

    for chain in unique_chains:

        this_chain = np.where(chains == chain)[0]
        other_chain = np.where(chains != chain)[0]

        for idx in this_chain:

            probs = contact_probs[idx, other_chain]

            score = float(np.max(probs))

            best_idx = other_chain[np.argmax(probs)]

            best_partner = f"{chains[best_idx]}{residues[best_idx]}"
            best_prob = float(np.max(probs))

            partners = []

            for j, p in zip(other_chain, probs):
                if p >= args.partner_threshold:
                    partners.append(
                        (
                            float(p),
                            f"{chains[j]}{residues[j]}({p:.2f})"
                        )
                    )

            # Sort partners from strongest to weakest
            partners.sort(reverse=True)

            partner_string = ";".join(
                partner for _, partner in partners
            )

            scores[(str(chain), int(residues[idx]))] = {
                "score": score,
                "best_partner": best_partner,
                "best_prob": best_prob,
                "n_partners": len(partners),
                "partners": partner_string
            }

    ############################################################
    # Read mmCIF
    ############################################################

    print("Reading mmCIF...")

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("AF3", args.cif)

    ############################################################
    # Apply scores to B-factors
    ############################################################

    print("Annotating structure...")

    output_rows = []

    for model in structure:

        for chain in model:

            for residue in chain:

                # Ignore hetero atoms/waters
                if residue.id[0] != " ":
                    continue

                chain_id = chain.id
                res_id = residue.id[1]

                entry = scores.get(
                    (chain_id, res_id),
                    {
                        "score": 0.0,
                        "best_partner": "",
                        "best_prob": 0.0,
                        "n_partners": 0,
                        "partners": ""
                    }
                )

                output_rows.append([
                    chain_id,
                    res_id,
                    residue.resname,
                    entry["score"],
                    entry["best_partner"],
                    entry["best_prob"],
                    entry["n_partners"],
                    entry["partners"]
                ])

                # Store interface score in B-factor (0-100)
                for atom in residue:
                    atom.set_bfactor(entry["score"] * 100)

    ############################################################
    # Write mmCIF
    ############################################################

    print("Writing annotated mmCIF...")

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(args.out_cif)

    ############################################################
    # Write TSV
    ############################################################

    print("Writing TSV...")

    with open(args.out_tsv, "w", newline="") as f:

        writer = csv.writer(f, delimiter="\t")

        writer.writerow([
            "chain",
            "residue",
            "aa",
            "interface_score",
            "best_partner",
            "best_prob",
            "n_partners",
            "partners"
        ])

        writer.writerows(output_rows)

    print("\nDone.")
    print(f"Annotated CIF : {args.out_cif}")
    print(f"Residue table : {args.out_tsv}")


if __name__ == "__main__":
    main()