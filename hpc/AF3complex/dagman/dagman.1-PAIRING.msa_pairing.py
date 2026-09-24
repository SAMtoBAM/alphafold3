#!/usr/bin/env python3
import os
import json
import tarfile
import argparse
import tempfile
import math
from pathlib import Path
import io


def load_probe(probe_path):
    """Load probe JSON from .json or .data_pipeline.tar.gz."""
    probe_path = Path(probe_path)

    if probe_path.name.endswith(".tar.gz"):
        with tempfile.TemporaryDirectory() as tmpdir:
            with tarfile.open(probe_path, "r:gz") as tar:
                members = tar.getmembers()
                json_member = next(
                    (m for m in members if m.name.endswith("_data.json")),
                    None
                )

                if json_member is None:
                    raise FileNotFoundError(
                        f"No *_data.json in {probe_path}"
                    )

                tar.extract(json_member, tmpdir)
                json_path = Path(tmpdir) / json_member.name

                with open(json_path) as f:
                    return json.load(f)

    else:
        with open(probe_path) as f:
            return json.load(f)


def combine_probe_msas(msa_dir, probe_path, out_dir, batch_size, highaccuracy):

    os.makedirs(out_dir, exist_ok=True)

    # ---------------------------------------------------------
    # Load probe once
    # ---------------------------------------------------------
    probe_json = load_probe(probe_path)
    probe_name = probe_json["name"]

    # ---------------------------------------------------------
    # Seeds
    # ---------------------------------------------------------
    if highaccuracy:
        seeds = [
            1, 11, 111, 1111, 11111,
            3, 33, 333, 3333, 33333,
            6, 66, 666, 6666, 66666,
            9, 99, 999, 9999, 99999
        ]
    else:
        seeds = [1]

    # ---------------------------------------------------------
    # Create job0 containing the probe protein by itself
    # ---------------------------------------------------------
    job0_dir = Path(out_dir) / "job0" / "inference_inputs"
    job0_dir.mkdir(parents=True, exist_ok=True)

    probe_only_tar = job0_dir / f"{probe_name}.data_pipeline.tar.gz"

    if not probe_only_tar.exists():

        # Copy probe protein
        probe_protein = probe_json["sequences"][0]["protein"].copy()

        # Force chain ID to A
        probe_protein["id"] = "A"

        # No paired MSA for standalone probe
        probe_protein["pairedMsa"] = ""

        # Build probe-only AF3 input
        probe_only = {
            "dialect": "alphafold3",
            "version": 1,
            "name": probe_name,
            "sequences": [
                {
                    "protein": probe_protein
                }
            ],
            "modelSeeds": seeds
        }

        # Write JSON directly into tar.gz
        data = json.dumps(
            probe_only,
            indent=2
        ).encode()

        tarinfo = tarfile.TarInfo(
            name=f"{probe_name}_data.json"
        )

        tarinfo.size = len(data)

        with tarfile.open(
            probe_only_tar,
            "w:gz"
        ) as out_tar:

            out_tar.addfile(
                tarinfo,
                io.BytesIO(data)
            )

        print(f"✔ job0: {probe_only_tar.name}")

    else:
        print(f"⏭ skip {probe_only_tar.name}")

    # ---------------------------------------------------------
    # Collect target inputs
    # ---------------------------------------------------------
    tarballs = sorted(
        Path(msa_dir).glob("*.data_pipeline.tar.gz")
    )

    if not tarballs:
        print(f"[!] No tarballs in {msa_dir}")
        print("[DONE] job0 created; no paired jobs generated")
        return

    total = len(tarballs)
    num_batches = math.ceil(total / batch_size)

    print(
        f"[INFO] {total} inputs → "
        f"{num_batches} paired job folders"
    )

    # ---------------------------------------------------------
    # Process paired batches
    # ---------------------------------------------------------
    for batch_idx in range(num_batches):

        start = batch_idx * batch_size
        end = min(start + batch_size, total)

        batch_files = tarballs[start:end]

        job_dir = (
            Path(out_dir)
            / f"job{batch_idx + 1}"
            / "inference_inputs"
        )

        job_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        # ONE temp dir per batch
        with tempfile.TemporaryDirectory() as tmpdir:

            for tar_path in batch_files:

                # Open target tar ONCE
                with tarfile.open(
                    tar_path,
                    "r:gz"
                ) as tar:

                    json_member = next(
                        (
                            m
                            for m in tar.getmembers()
                            if m.name.endswith("_data.json")
                        ),
                        None
                    )

                    if json_member is None:
                        print(
                            f"[WARN] missing JSON in {tar_path}"
                        )
                        continue

                    # Extract only required JSON
                    tar.extract(
                        json_member,
                        tmpdir
                    )

                    extracted_json = (
                        Path(tmpdir) / json_member.name
                    )

                # Load target JSON
                with open(extracted_json) as f:
                    target_json = json.load(f)

                target_name = target_json["name"]

                expected_tar = (
                    job_dir
                    / f"{probe_name}_{target_name}"
                      f".data_pipeline.tar.gz"
                )

                if expected_tar.exists():
                    print(
                        f"⏭ skip {expected_tar.name}"
                    )
                    continue

                # Build combined structure
                combined = {
                    "dialect": "alphafold3",
                    "version": 1,
                    "name": f"{probe_name}_{target_name}",
                    "sequences": [],
                    "modelSeeds": seeds
                }

                # Probe chain
                probe_protein = (
                    probe_json["sequences"][0]["protein"].copy()
                )

                probe_protein["id"] = "A"
                probe_protein["pairedMsa"] = ""

                combined["sequences"].append(
                    {
                        "protein": probe_protein
                    }
                )

                # Target chain
                target_protein = (
                    target_json["sequences"][0]["protein"].copy()
                )

                target_protein["id"] = "B"
                target_protein["pairedMsa"] = ""

                combined["sequences"].append(
                    {
                        "protein": target_protein
                    }
                )

                combined_json_name = (
                    f"{probe_name}_{target_name}_data.json"
                )

                # Write directly into tar.gz
                data = json.dumps(
                    combined,
                    indent=2
                ).encode()

                tarinfo = tarfile.TarInfo(
                    name=combined_json_name
                )

                tarinfo.size = len(data)

                with tarfile.open(
                    expected_tar,
                    "w:gz"
                ) as out_tar:

                    out_tar.addfile(
                        tarinfo,
                        io.BytesIO(data)
                    )

                print(
                    f"✔ job{batch_idx + 1}: "
                    f"{expected_tar.name}"
                )

    print(
        f"[DONE] {total} pairs across "
        f"{num_batches} jobs + job0 probe-only"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--msa_dir",
        required=True
    )

    parser.add_argument(
        "--probe",
        required=True
    )

    parser.add_argument(
        "--out_dir",
        required=True
    )

    parser.add_argument(
        "--batch_size",
        type=int,
        default=10
    )

    parser.add_argument(
        "--highaccuracy",
        action="store_true"
    )

    args = parser.parse_args()

    combine_probe_msas(
        args.msa_dir,
        args.probe,
        args.out_dir,
        args.batch_size,
        args.highaccuracy
    )