#!/usr/bin/env python3
import os
import json
import tarfile
import argparse
import tempfile
import math
from pathlib import Path


def load_probe(probe_path):
    """Load probe JSON from either plain .json or .data_pipeline.tar.gz."""
    probe_path = Path(probe_path)
    if probe_path.suffixes[-2:] == [".tar", ".gz"] or probe_path.name.endswith(".data_pipeline.tar.gz"):
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


def combine_probe_msas(msa_dir, probe_path, out_dir, batch_size):
    os.makedirs(out_dir, exist_ok=True)

    # Load probe JSON once
    probe_json = load_probe(probe_path)
    probe_name = probe_json["name"]

    # Find input tarballs
    tarballs = sorted(Path(msa_dir).glob("*.data_pipeline.tar.gz"))
    if not tarballs:
        print(f"[!] No .data_pipeline.tar.gz files found in {msa_dir}")
        return

    total = len(tarballs)
    num_batches = math.ceil(total / batch_size)
    print(f"Found {total} inputs → {num_batches} job folders (batch size = {batch_size})")

    # Process in batches
    for batch_idx in range(num_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, total)
        batch_files = tarballs[start:end]

        job_dir = Path(out_dir) / f"job{batch_idx + 1}" / "inference_inputs"
        job_dir.mkdir(parents=True, exist_ok=True)

        for tar_path in batch_files:
            with tempfile.TemporaryDirectory() as tmpdir:
                # Extract input tarball
                with tarfile.open(tar_path, "r:gz") as tar:
                    tar.extractall(tmpdir)

                # Find contained *_data.json
                json_files = list(Path(tmpdir).glob("*_data.json"))
                if not json_files:
                    print(f"[!] No *_data.json inside {tar_path}")
                    continue

                target_json_path = json_files[0]
                with open(target_json_path) as f:
                    target_json = json.load(f)

                target_name = target_json["name"]

                # Build combined JSON with modelSeeds = [1]
                combined = {
                    "dialect": "alphafold3",
                    "version": 1,
                    "name": f"{probe_name}_{target_name}",
                    "sequences": [],
                    "modelSeeds": [1]  # <-- required for AlphaFold3
                }

                # Chain A → probe
                probe_protein = probe_json["sequences"][0]["protein"].copy()
                probe_protein["id"] = "A"
                probe_protein["pairedMsa"] = ""
                combined["sequences"].append({"protein": probe_protein})

                # Chain B → target
                target_protein = target_json["sequences"][0]["protein"].copy()
                target_protein["id"] = "B"
                target_protein["pairedMsa"] = ""
                combined["sequences"].append({"protein": target_protein})

                # Write combined JSON
                combined_json_name = f"{probe_name}_{target_name}_data.json"
                combined_json_path = job_dir / combined_json_name
                with open(combined_json_path, "w") as f:
                    json.dump(combined, f, indent=2)

                # Tarball it
                tar_output = job_dir / f"{probe_name}_{target_name}.data_pipeline.tar.gz"
                with tarfile.open(tar_output, "w:gz") as tar:
                    tar.add(combined_json_path, arcname=combined_json_name)

                os.remove(combined_json_path)
                print(f"✔ job{batch_idx + 1}: {tar_output.name}")

    print(f"Finished generating {total} paired MSAs across {num_batches} job folders.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine single-protein AlphaFold3 MSA JSONs with a probe JSON (or tarball) into batched job directories."
    )
    parser.add_argument("--msa_dir", required=True, help="Directory with *.data_pipeline.tar.gz MSA archives")
    parser.add_argument("--probe", required=True, help="Path to probe JSON or .data_pipeline.tar.gz")
    parser.add_argument("--out_dir", required=True, help="Base output directory for job folders")
    parser.add_argument("--batch_size", type=int, default=10, help="Number of tarballs per job folder (default: 10)")

    args = parser.parse_args()
    combine_probe_msas(args.msa_dir, args.probe, args.out_dir, args.batch_size)

