#!/usr/bin/env python3
import os
import json
import tarfile
import argparse
import math
from pathlib import Path
from Bio import SeqIO  # pip install biopython


def load_probe(probe_path):
    """Load probe JSON from either plain .json or .data_pipeline.tar.gz."""
    probe_path = Path(probe_path)
    if probe_path.suffixes[-2:] == [".tar", ".gz"] or probe_path.name.endswith(".data_pipeline.tar.gz"):
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            with tarfile.open(probe_path, "r:gz") as tar:
                tar.extractall(tmpdir)
            json_files = list(Path(tmpdir).glob("*_data.json"))
            if not json_files:
                raise FileNotFoundError(f"No *_data.json found inside {probe_path}")
            with open(json_files[0]) as f:
                probe_json = json.load(f)
        return probe_json
    else:
        with open(probe_path) as f:
            return json.load(f)


def combine_probe_targets(probe_path, fasta_path, out_dir, batch_size):
    os.makedirs(out_dir, exist_ok=True)

    # Load probe JSON
    probe_json = load_probe(probe_path)
    probe_name = probe_json["name"]

    # Load nucleotide targets
    fasta_path = Path(fasta_path)
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"[!] No sequences found in {fasta_path}")
        return

    total = len(records)
    num_batches = math.ceil(total / batch_size)
    print(f"Found {total} sequences → {num_batches} job folders (batch size = {batch_size})")

    for batch_idx in range(num_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, total)
        batch_records = records[start:end]

        job_dir = Path(out_dir) / f"job{batch_idx + 1}" / "inference_inputs"
        job_dir.mkdir(parents=True, exist_ok=True)

        for record in batch_records:
            target_name = record.id
            target_seq = str(record.seq)

            # Build combined JSON
            combined = {
                "dialect": "alphafold3",
                "version": 1,
                "name": f"{probe_name}_{target_name}",
                "sequences": [],
                "modelSeeds": [1],
                "bondedAtomPairs": None,
                "userCCD": None
            }

            # Chain A → probe (protein, fully preserved)
            probe_protein = json.loads(json.dumps(probe_json["sequences"][0]["protein"]))  # deep copy
            probe_protein["id"] = "A"
            probe_protein["pairedMsa"] = ""
            combined["sequences"].append({"protein": probe_protein})

            # Chain B → DNA (new)
            target_dna = {
                "id": "B",
                "sequence": target_seq,
                "modifications": []
            }
            combined["sequences"].append({"dna": target_dna})

            # Write combined JSON
            json_name = f"{probe_name}_{target_name}_data.json"
            json_path = job_dir / json_name
            with open(json_path, "w") as f:
                json.dump(combined, f, indent=2)

            # Tarball the JSON
            tar_name = f"{probe_name}_{target_name}.data_pipeline.tar.gz"
            tar_path = job_dir / tar_name
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(json_path, arcname=json_name)

            os.remove(json_path)
            print(f"✔ job{batch_idx + 1}: {tar_name}")

    print(f"Finished generating {total} paired sequences across {num_batches} job folders.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine a probe JSON with each DNA sequence from a multi-FASTA file into AlphaFold3 job JSONs (no MSA), preserving template/query indices."
    )
    parser.add_argument("--probe", required=True, help="Path to probe JSON or .data_pipeline.tar.gz")
    parser.add_argument("--nucleotides", required=True, help="Multi-FASTA file with DNA sequences")
    parser.add_argument("--out_dir", required=True, help="Output directory for JSON/tarball files")
    parser.add_argument("--batch_size", type=int, default=10, help="Number of sequences per job folder (default: 10)")

    args = parser.parse_args()
    combine_probe_targets(args.probe, args.nucleotides, args.out_dir, args.batch_size)

