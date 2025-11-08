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


def load_fasta_records(fasta_path):
    """Load sequences from a FASTA file."""
    fasta_path = Path(fasta_path)
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No entries found in {fasta_path}")
    return records


def combine_probe_with_targets(probe_path, fasta_path, out_dir, batch_size, mode):
    os.makedirs(out_dir, exist_ok=True)

    # Load probe JSON
    probe_json = load_probe(probe_path)
    probe_name = probe_json["name"]

    # Load FASTA entries
    records = load_fasta_records(fasta_path)
    total = len(records)
    num_batches = math.ceil(total / batch_size)
    print(f"Found {total} {mode} entries → {num_batches} job folders (batch size = {batch_size})")

    for batch_idx in range(num_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, total)
        batch_records = records[start:end]

        job_dir = Path(out_dir) / f"job{batch_idx + 1}" / "inference_inputs"
        job_dir.mkdir(parents=True, exist_ok=True)

        for record in batch_records:
            target_name = record.id
            target_seq = str(record.seq)

            combined = {
                "dialect": "alphafold3",
                "version": 1,
                "name": f"{probe_name}_{target_name}",
                "sequences": [],
                "modelSeeds": [1],
                "bondedAtomPairs": None,
                "userCCD": None
            }

            # Chain A → protein (probe)
            probe_protein = json.loads(json.dumps(probe_json["sequences"][0]["protein"]))  # deep copy
            probe_protein["id"] = "A"
            probe_protein["pairedMsa"] = ""
            combined["sequences"].append({"protein": probe_protein})

            # Chain B → either DNA or ligand
            if mode == "nucleotides":
                combined["sequences"].append({
                    "dna": {
                        "id": "B",
                        "sequence": target_seq,
                        "modifications": []
                    }
                })
            elif mode == "ligand":
                combined["sequences"].append({
                    "ligand": {
                        "id": "B",
                        "smiles": target_seq
                    }
                })
            else:
                raise ValueError("Invalid mode — must be 'nucleotides' or 'ligand'")

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

    print(f"Finished generating {total} {mode} combinations across {num_batches} job folders.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine a protein probe with either DNA sequences OR ligands (mutually exclusive) into AlphaFold3 job JSONs."
    )
    parser.add_argument("--probe", required=True, help="Path to probe JSON or .data_pipeline.tar.gz")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--nucleotides", help="Multi-FASTA file with DNA sequences")
    group.add_argument("--ligands", help="FASTA file with ligands (SMILES strings)")
    parser.add_argument("--out_dir", required=True, help="Output directory for JSON/tarball files")
    parser.add_argument("--batch_size", type=int, default=10, help="Number of sequences per job folder (default: 10)")

    args = parser.parse_args()

    if args.nucleotides:
        mode = "nucleotides"
        fasta_path = args.nucleotides
    else:
        mode = "ligand"
        fasta_path = args.ligands

    combine_probe_with_targets(args.probe, fasta_path, args.out_dir, args.batch_size, mode)
