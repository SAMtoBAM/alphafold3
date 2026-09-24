#!/usr/bin/env bash
set -euo pipefail

probe="${1,,}"   # lowercase safeguard

tar -xzf proteome_msa.tar.gz -C ./

python dagman.1-PAIRING.msa_pairing.py \
  --msa_dir ./proteome_msa/ \
  --probe "./proteome_msa/${probe}.data_pipeline.tar.gz" \
  --out_dir ./ \
  --batch_size 20

tar -czf msa_pairing_output.tar.gz job*
