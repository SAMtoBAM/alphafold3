#!/usr/bin/env python3
"""
Given a probe FASTA and a proteome FASTA, this script will create job folders (job1, job2, …)
each containing AlphaFold Server JSON files (AF3Complex-compatible) describing complexes of:
    probe + each protein in the proteome

Usage:
    python set_up_directory_complex.py probe.fasta proteome.fasta batch_size

Example:
    python set_up_directory_complex.py probe.fasta proteome.fasta 10
"""

import json
from Bio import SeqIO
import os
import argparse
import math


def combine_probe_and_proteome(probe_fasta, proteome_fasta, batch_size):
    # Read probe (expecting 1 sequence)
    probe_records = list(SeqIO.parse(probe_fasta, "fasta"))
    if len(probe_records) != 1:
        raise ValueError("Probe FASTA must contain exactly one sequence.")
    probe_record = probe_records[0]
    probe_seq = str(probe_record.seq).strip()

    # Read proteome sequences
    proteome_records = list(SeqIO.parse(proteome_fasta, "fasta"))
    total = len(proteome_records)
    num_batches = math.ceil(total / batch_size)

    print(f"Loaded 1 probe ({probe_record.id}) and {total} proteome sequences.")
    print(f"→ Creating {num_batches} job folders with up to {batch_size} complexes each.")

    def build_complex_entry(proteome_record):
        """Build one AlphaFold3 complex entry (probe + proteome sequence)."""
        prot_seq = str(proteome_record.seq).strip()

        seq_list = [
            {"proteinChain": {"sequence": probe_seq, "count": 1}},     # probe
            {"proteinChain": {"sequence": prot_seq, "count": 1}}       # target protein
        ]

        entry = {
            "name": f"{probe_record.id}_{proteome_record.id}",
            "modelSeeds": [1, 11, 111, 1111, 11111, 3, 33, 333, 3333, 33333, 6, 66, 666, 6666, 66666, 9, 99, 999, 9999, 99999],
            "sequences": seq_list,
            "dialect": "alphafoldserver",
            "version": 1
        }
        return entry

    # Create job folders and write JSONs
    for batch_idx in range(num_batches):
        job_dir = f"job{batch_idx + 1}"
        data_inputs_dir = os.path.join(job_dir, "data_inputs")
        inference_inputs_dir = os.path.join(job_dir, "inference_inputs")
        os.makedirs(data_inputs_dir, exist_ok=True)
        os.makedirs(inference_inputs_dir, exist_ok=True)

        # Select batch
        batch_records = proteome_records[batch_idx * batch_size:(batch_idx + 1) * batch_size]
        batch_entries = [build_complex_entry(r) for r in batch_records]

        # Write combined JSON
        json_path = os.path.join(data_inputs_dir, "fold_input_server.json")
        with open(json_path, "w") as f:
            json.dump(batch_entries, f, indent=2)

        print(f"✔ Wrote {len(batch_records)} complexes to {json_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine a probe and proteome FASTA into batched AlphaFold Server JSONs for complex prediction."
    )
    parser.add_argument("probe_fasta", type=str, help="Path to probe FASTA (one sequence).")
    parser.add_argument("proteome_fasta", type=str, help="Path to proteome FASTA (many sequences).")
    parser.add_argument("batch_size", type=int, help="Number of complexes per job folder.")

    args = parser.parse_args()
    combine_probe_and_proteome(args.probe_fasta, args.proteome_fasta, args.batch_size)
