#!/usr/bin/env bash

exec > inference_check.debug.log 2>&1

set -euo pipefail

shopt -s nullglob


outputs=(job*/*.inference_pipeline.tar.gz)
inputs=(job*/inference_inputs/*.data_pipeline.tar.gz)

outputcount=${#outputs[@]}
inputcount=${#inputs[@]}

if [[ "$outputcount" -eq "$inputcount" ]]; then
    echo "Looks good! Number of output files match input (${inputcount})"
else
    echo "Number of input inference files does not match number of output folders"
    echo "Output count: ${outputcount}"
    echo "Input count: ${inputcount}"
    exit 1
fi

mkdir proteome_inference
mv job*/*inference_pipeline.tar.gz proteome_inference/


mkdir rerun-job
mkdir rerun-job/inference_inputs
missing_cif="0"
for f in proteome_inference/*.inference_pipeline.tar.gz; do
    [[ -e "$f" ]] || continue

    if ! grep -iq '\.cif$' < <(tar -tf "$f" 2>/dev/null); then
        echo "WARNING: no .cif found in $f"
        ((++missing_cif))

        complex=$(basename "$f" .inference_pipeline.tar.gz)

        src=$(find job*/inference_inputs \
            -maxdepth 1 \
            -type f \
            -iname "${complex}.data_pipeline.tar.gz" \
            -print -quit)

        if [[ -z $src ]]; then
            echo "WARNING: no data_pipeline tar found for $complex"
            continue
        fi

        cp "$src" rerun-job/inference_inputs/
    fi
done
[[ "$missing_cif" -eq 0 ]] && echo "All inference_pipeline archives contain .cif files"
if [ "$missing_cif" != 0 ]
then
echo "Empty results folders are usually due to an issue with VRAM available by the GPU"
echo "Therefore creating a new sub file = 'inference_pipeline.complex.rerun.sub'"
echo "The new file has a higher VRAM minimum (from 48GB to 80GB) and only runs on only the data in the folder 'job-rerun'"
cat inference_pipeline.complex.sub | sed 's/CUDAGlobalMemoryMb > 48000/CUDAGlobalMemoryMb > 80000/' | sed 's/job\*/rerun-job/' > inference_pipeline.complex.rerun.sub
fi

##resubmit the failed inferences with the higher VRAM minimum request
#condor_submit inference_pipeline.complex.rerun.sub
##now move the output into the proteome_inference folder again (overwriting the old empty ones)
#mv rerun-job/*inference_pipeline.tar.gz proteome_inference/
#rm -r rerun-job/

##now run the check above again!
##if this didn't work, likely you will need to split the protein (trying not to break domains); re-run the MSA step; then rerun inference. or just skip it.


##once finished, tidy up the input data
#rm -r job*

