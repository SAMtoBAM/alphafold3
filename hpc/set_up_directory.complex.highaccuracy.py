# Given a multifasta, this script will create a folder named job1...jobN with subfolders named data_inputs and inference_inputs.
# This script will create a custom json file for each fasta sequence, and place them in batches of batch size (e.g. 10) under each job1 to jobN folders
# the batch size can be provided by the user through a command line argument.
# usage: set_up_folders.py [-h] fasta_file batch_size

# for example, let test.fasta look like this, and batch_size = 5:
# >sequence1
# ABCDEFG
# >sequence2
# HIJKLMNOP
# >sequence3
# QRSTUV
# >sequence4
# ABCDEFG
# >sequence5
# HIJKLMNOP
# >sequence6
# QRSTUV
# This script will create these 2 job folders, with 5 json files in each.
# job1
# ├── data_inputs
# │   └── fold_input_1.json
# │   └── fold_input_2.json
# │   └── fold_input_3.json
# │   └── fold_input_4.json
# │   └── fold_input_5.json
# └── inference_inputs
# job2
# ├── data_inputs
# │   └── fold_input_6.json
# └── inference_inputs

## USAGE ##
# chmod +x set_up_directory.py
# python set_up_directory test.fasta batch_size

# Import libraries

import json
from Bio import SeqIO
import os
import argparse
import math


def fasta_to_server_json(fasta_file, batch_size):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    count = len(sequences)
    num_batches = math.ceil(count / batch_size)

    def build_entry(record):
        """Build one entry for AlphaFold Server format (compatible with AF3Complex)."""
        seq_str = str(record.seq).strip()

        # Detect complex by ':' separator in sequence
        if ':' in seq_str:
            sub_sequences = seq_str.split(':')
            seq_list = [
                {
                    "proteinChain": {
                        "sequence": subseq,
                        "count": 1
                    }
                }
                for subseq in sub_sequences
            ]
            print(f"[Complex] {record.id} → {len(sub_sequences)} chains")
        else:
            seq_list = [
                {
                    "proteinChain": {
                        "sequence": seq_str,
                        "count": 1
                    }
                }
            ]
            print(f"[Single]  {record.id}")

        # Build entry following AlphaFold Server JSON dialect
        entry = {
            "name": record.id,
            "modelSeeds": [1, 11, 111, 1111, 11111],
            "sequences": seq_list,
            "dialect": "alphafoldserver",
            "version": 1
        }
        return entry

    # Create job folders
    for batch_index in range(num_batches):
        job_dir = f"job{batch_index + 1}"
        data_inputs_dir = os.path.join(job_dir, "data_inputs")
        inference_inputs_dir = os.path.join(job_dir, "inference_inputs")
        os.makedirs(data_inputs_dir, exist_ok=True)
        os.makedirs(inference_inputs_dir, exist_ok=True)

        # Records for this job
        batch_records = sequences[batch_index * batch_size:(batch_index + 1) * batch_size]
        batch_json_list = [build_entry(r) for r in batch_records]

        # Write one JSON file per batch (AlphaFold Server format → list at top level)
        json_filename = os.path.join(data_inputs_dir, "fold_input_server.json")
        with open(json_filename, "w") as json_file:
            json.dump(batch_json_list, json_file, indent=2)

        print(f"✔ Wrote {len(batch_records)} entries → {json_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA to per-job AlphaFold Server JSON format (AF3Complex compatible).")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file")
    parser.add_argument("batch_size", type=int, help="Number of sequences per job folder")

    args = parser.parse_args()
    fasta_to_server_json(args.fasta_file, args.batch_size)

