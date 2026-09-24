#!/usr/bin/env bash
set -euo pipefail

shopt -s nullglob

source ~/miniconda3/etc/profile.d/conda.sh
conda activate AF3

files=(rerun-job/*inference_pipeline.tar.gz)

if (( ${#files[@]} )); then
    mv "${files[@]}" proteome_inference/
fi

##clean up previous files
rm -rf job* rerun-job

##remove debug log if succeeded with inference steps
[ -f inference_check.debug.log ] && rm inference_check.debug.log

##compress output for next step as it needs to be transferred to the submission node again
tar -czf proteome_inference.tar.gz proteome_inference &&
rm -r proteome_inference
