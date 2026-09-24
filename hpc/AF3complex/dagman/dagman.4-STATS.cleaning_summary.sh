#!/usr/bin/env bash
set -euo pipefail

probe="${1,,}"   # lowercase safeguard

##unzip and activate the conda envs
tar -xzf AF3.conda.tar.gz
export PATH="$PWD/bin:$PATH"

echo "probe: ${probe}"

echo "unpacking inference folder"
tar -xzf proteome_inference.tar.gz -C ./ && rm proteome_inference.tar.gz


#echo "cleaning up alphafold3 output to keep only necessary files"
#mkdir proteome_inference_clean
#for f in proteome_inference/*.inference_pipeline.tar.gz
#do
#  complex=$(basename "${f}" .inference_pipeline.tar.gz)
#  #echo "Processing $complex"
#  tmpdir=$(mktemp -d)
#  # Extract only kept files into temp dir
#  tar -xzf "${f}" --exclude='seed*' --exclude='*_data.json' -C "$tmpdir"
#  # Repack
#  tar -czf "proteome_inference_clean/${complex}.inference_pipeline.tar.gz" -C "$tmpdir" .
#  rm -rf "$tmpdir"
#done

#rm -r proteome_inference
#mv proteome_inference_clean proteome_inference

echo "getting per residue pLDDT and ipSAE values"

##header for the ipSAE values
echo "probe;pair;ipSAE;ipSAE_d0chn;ipSAE_d0dom" | tr ';' '\t' > "${probe}_pairs.ipSAE.tsv"


##get the pLDDT values for the single probe by itself
mkdir ${probe}.pLDDT

tar -xzf proteome_inference/${probe}.inference_pipeline.tar.gz \
        "./${probe}_confidences.json" \
        "./${probe}_model.cif" \
        "./${probe}_summary_confidences.json"

python3 dagman.4-STATS.extract_probe_plddt.py \
    "${probe}_confidences.json" \
    "${probe}_model.cif" \
    > ${probe}.pLDDT/${probe}.residue_plddt.tsv

rm ${probe}_summary_confidences.json ${probe}_model.cif ${probe}_confidences.json



##now loop through each pairing and calculate the ipSAE metrics and put them in the temp tsv
while read -r pair
do
    folder=$(find proteome_inference \
        -maxdepth 1 \
        -type f \
        -iname "${probe}_${pair}.inference_pipeline.tar.gz" \
        -print -quit)

    if [[ -z "$folder" ]]; then
        echo "WARNING: no tar.gz found for ${probe}_${pair}"
        echo "${probe};${pair};NA;NA;NA" | tr ';' '\t' >> "${probe}_pairs.ipSAE.tsv"
        continue
    fi

    protein=$(basename "$folder" .inference_pipeline.tar.gz)

    echo "Processing ${protein}"

    tar -xzf "$folder" \
        "./${protein}_confidences.json" \
        "./${protein}_model.cif" \
        "./${protein}_summary_confidences.json" \
        2>/dev/null || true

    if [[ ! -f "${protein}_model.cif" ]]; then
        echo "WARNING: no CIF found for ${protein}"
        echo "${probe};${pair};NA;NA;NA" | tr ';' '\t' >> "${probe}_pairs.ipSAE.tsv"

        rm -f "${protein}_confidences.json" \
              "${protein}_model.cif" \
              "${protein}_summary_confidences.json"

        continue
    fi

    af3tools metrics \
        -p "${protein}_confidences.json" \
        -s "${protein}_model.cif"

    awk -v probe="$probe" -v pair="$pair" \
        '$5 == "max" {
            print probe "\t" pair "\t" $6 "\t" $7 "\t" $8
        }' \
        "${protein}_model_10_10.txt" \
        >> "${probe}_pairs.ipSAE.tsv"



    python3 dagman.4-STATS.extract_probe_plddt.py \
    "${protein}_confidences.json" \
    "${protein}_model.cif" \
    > ${probe}.pLDDT/${protein}.residue_plddt.tsv

    rm -f "${protein}_confidences.json" \
          "${protein}_model.cif" \
          "${protein}_summary_confidences.json" \
          "${protein}_model_10_10"*

done < full_protein_list.txt


##combine all the files into a single output
##first get a list of the files with the single probe being explicitely the first
files=(
    "${probe}.pLDDT/${probe}.residue_plddt.tsv"
)

for file in "${probe}.pLDDT/${probe}_"*.residue_plddt.tsv; do
    files+=("$file")
done

##now combine all the files, grabbing all the columns for the probe alone
##then only the third column for the rest
awk '
BEGIN {
    OFS="\t"
}

FNR == 1 {
    file[++n] = FILENAME
    next
}

{
    residue = $1
    aa[residue] = $2
    plddt[FILENAME,residue] = $3

    if (!(residue in seen)) {
        residues[++nr] = residue
        seen[residue] = 1
    }
}

END {

    printf "residue\taa"

    for (i = 1; i <= n; i++) {
        name = file[i]
        sub(/^.*\//, "", name)
        sub(/\.residue_plddt\.tsv$/, "", name)

        printf "%s%s", OFS, name
    }

    printf "\n"

    for (r = 1; r <= nr; r++) {

        residue = residues[r]

        printf "%s%s%s", residue, OFS, aa[residue]

        for (i = 1; i <= n; i++) {
            printf "%s%s", OFS, plddt[file[i],residue]
        }

        printf "\n"
    }
}
' "${files[@]}" > "${probe}.combined_plddt.tsv"


##Assuming everything ran well we now want to evaluate the best model-sample combination per comple
##already the best seed has been evaluated for the last model/seed run. So if only one seed was given (default) then run the below getting information on the best seed for that seed
##this handles getting confidence scores whether or not the file is completed and fills in NAs if nothing
##folder search also needs to allow for some letter being made lowercase by alphafold

echo "Generating a summary file of the scores per probe-pair complex"

echo "protein1;protein2;fraction_disordered;has_clash;pLDDT;ptm;iptm;ranking_score" | tr ';' '\t' > confidence_summary.tsv

##old version with jq used
#cat full_protein_list.txt | while read -r pair
#do
#    folder=$(find proteome_inference \
#    -maxdepth 1 \
#    -type f \
#    -iname "${probe}_${pair}.inference_pipeline.tar.gz" \
#    -print -quit)
#    protein=$(basename "$folder" .inference_pipeline.tar.gz)
#
#    ## defaults
#    pLDDT="NA"
#    iptm="NA"
#    ptm="NA"
#    RS="NA"
#    FD="NA"
#    HC="NA"
#
#    ## pLDDT
#    cif_out=$(tar -axf "$folder" "./${protein}_model.cif" -O 2>/dev/null)
#    if [[ -n "$cif_out" ]]; then
#        pLDDT=$(echo "$cif_out" \
#            | awk '/_ma_qa_metric_global.metric_value/ {print $2; exit}')
#    fi
#
#    ## JSON metrics (SAFE)
#    if tar -tf "$folder" "./${protein}_summary_confidences.json" &>/dev/null; then
#        read iptm ptm RS FD HC < <(
#            tar -axf "$folder" "./${protein}_summary_confidences.json" -O \
#            | jq -r '
#                [
#                  .iptm,
#                  .ptm,
#                  .ranking_score,
#                  .fraction_disordered,
#                  .has_clash
#                ] | map(. // "NA") | @tsv
#            '
#        )
#    fi
#
#    ## output
#    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
#        "$probe" "$pair" "$FD" "$HC" "$pLDDT" "$ptm" "$iptm" "$RS"
#
#done >> confidence_summary.tsv

echo "Generating a summary file of the scores per probe-pair complex again"

ls

while read -r pair
do
    folder=$(find proteome_inference \
        -maxdepth 1 \
        -type f \
        -iname "${probe}_${pair}.inference_pipeline.tar.gz" \
        -print -quit)

    ## defaults
    pLDDT="NA"
    iptm="NA"
    ptm="NA"
    RS="NA"
    FD="NA"
    HC="NA"

    ## only extract if archive exists
    if [[ -n "$folder" ]]; then

        protein=$(basename "$folder" .inference_pipeline.tar.gz)

        ## pLDDT
        pLDDT=$(
            tar -axf "$folder" "./${protein}_model.cif" -O 2>/dev/null \
            | awk '/_ma_qa_metric_global.metric_value/ {print $2}' \
            | head -n1
        ) || pLDDT="NA"

        [[ -z "$pLDDT" ]] && pLDDT="NA"

        ## JSON metrics
        if tar -tf "$folder" "./${protein}_summary_confidences.json" &>/dev/null; then

            read iptm ptm RS FD HC < <(
                tar -axf "$folder" "./${protein}_summary_confidences.json" -O \
                | python3 -c '
import json
import sys

try:
    d=json.load(sys.stdin)
    print(
        d.get("iptm","NA"),
        d.get("ptm","NA"),
        d.get("ranking_score","NA"),
        d.get("fraction_disordered","NA"),
        d.get("has_clash","NA"),
        sep="\t"
    )
except Exception:
    print("NA\tNA\tNA\tNA\tNA")
'
            )
        fi
    fi

    ## always output one row
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$probe" \
        "$pair" \
        "$FD" \
        "$HC" \
        "$pLDDT" \
        "$ptm" \
        "$iptm" \
        "$RS"

done < full_protein_list.txt >> confidence_summary.tsv


##combine the ipSAE values with all the other metrics
cut -f3-5 ${probe}_pairs.ipSAE.tsv | paste confidence_summary.tsv - > temp
mv temp confidence_summary.tsv

##now we can use these confidence scores to try and remove some really unlikely complexes in order to save room in storage
##saving all proteins = ~10GB per probe...therefore we would have 10k*10GB=100TB
##so if we remove ~99% of bad calls we can reduce this to around 1TB...more reasonable 

##use emperical p-values to find the good candidates (also use the ranking score and the iptm seperately and take the nonredundant list of candidates)
##with a cut of of 0.01 (lenient here in order to be forgiving but save ALOT of space)

echo "Getting the top 1% of scores and only keeping models for these complexes"

##Step 1: Count total number of scores
total=$( tail -n+2 confidence_summary.tsv | wc -l)

##Step 2: Compute empirical p-values using the ranking score and grab those with a pval<0.05 (and ignoring NAs)
tail -n+2 confidence_summary.tsv |
awk -F"\t" '
$8 != "NA" {
    scores[++m] = $8
}
{
    lines[NR] = $0
    vals[NR] = $8
}
END {
    n = m
    for (i = 1; i <= NR; i++) {
        if (vals[i] == "NA") {
            print lines[i] "\tNA\tNA"
            continue
        }
        count = 0
        for (j = 1; j <= n; j++) {
            if (scores[j] <= vals[i]) count++
        }
        p = 1 - (count / n)
        label = (p < 0.01 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2 > temp_list.txt

##Step 3: same as above but using only the iptm
tail -n+2 confidence_summary.tsv |
awk -F"\t" '
$7 != "NA" {
    scores[++m] = $7
}
{
    lines[NR] = $0
    vals[NR] = $7
}
END {
    n = m
    for (i = 1; i <= NR; i++) {
        if (vals[i] == "NA") {
            print lines[i] "\tNA\tNA"
            continue
        }
        count = 0
        for (j = 1; j <= n; j++) {
            if (scores[j] <= vals[i]) count++
        }
        p = 1 - (count / n)
        label = (p < 0.01 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 4: same as above but using only the ipSAE
tail -n+2 confidence_summary.tsv |
awk -F"\t" '
$7 != "NA" {
    scores[++m] = $9
}
{
    lines[NR] = $0
    vals[NR] = $9
}
END {
    n = m
    for (i = 1; i <= NR; i++) {
        if (vals[i] == "NA") {
            print lines[i] "\tNA\tNA"
            continue
        }
        count = 0
        for (j = 1; j <= n; j++) {
            if (scores[j] <= vals[i]) count++
        }
        p = 1 - (count / n)
        label = (p < 0.01 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 5: same as above but using only the ipSAE_d0chn
tail -n+2 confidence_summary.tsv |
awk -F"\t" '
$7 != "NA" {
    scores[++m] = $10
}
{
    lines[NR] = $0
    vals[NR] = $10
}
END {
    n = m
    for (i = 1; i <= NR; i++) {
        if (vals[i] == "NA") {
            print lines[i] "\tNA\tNA"
            continue
        }
        count = 0
        for (j = 1; j <= n; j++) {
            if (scores[j] <= vals[i]) count++
        }
        p = 1 - (count / n)
        label = (p < 0.01 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 6: same as above but using only the ipSAE_d0dom
tail -n+2 confidence_summary.tsv |
awk -F"\t" '
$7 != "NA" {
    scores[++m] = $11
}
{
    lines[NR] = $0
    vals[NR] = $11
}
END {
    n = m
    for (i = 1; i <= NR; i++) {
        if (vals[i] == "NA") {
            print lines[i] "\tNA\tNA"
            continue
        }
        count = 0
        for (j = 1; j <= n; j++) {
            if (scores[j] <= vals[i]) count++
        }
        p = 1 - (count / n)
        label = (p < 0.01 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 7: get a nonredundant list but only keeping those that were candidates in at least three of the lists. i.e. found in the top of 3 different metrics.
awk '{count[$1]++; values[NR]=$1} END {for (i=1;i<=NR;i++) if (count[values[i]] >= 3) print values[i]}' temp_list.txt | sort -u > candidate_list.txt
rm temp_list.txt

##Step 5: now only keep data for those in the candidate_list.txt file
mkdir proteome_inference_candidates
cat candidate_list.txt | while read -r pair
do
    folder=$(find proteome_inference \
    -maxdepth 1 \
    -type f \
    -iname "${probe}_${pair}.inference_pipeline.tar.gz" \
    -print -quit)
    mv ${folder} proteome_inference_candidates/${probe}_${pair}.inference_pipeline.tar.gz
done

##also transfer probe by itself
mv proteome_inference/${probe}.inference_pipeline.tar.gz proteome_inference_candidates/${probe}.inference_pipeline.tar.gz

##get rid of the rest and rename
rm -r proteome_inference
mv proteome_inference_candidates proteome_inference

echo "Packing up the results"
##reduce the file numbers for storage by zipping up all the already zipped results
tar -czf proteome_inference.tar.gz proteome_inference && rm -r proteome_inference

##do after moving back to cluster
#cp confidence_summary.tsv confidence_summaries/${probe}.confidence_summary.tsv

cp confidence_summary.tsv ${probe}.confidence_summary.tsv


##pack up the output for the probe and move it off for temp storage
mkdir ${probe}_complexes
mv proteome_inference.tar.gz ${probe}_complexes
cp ${probe}.confidence_summary.tsv ${probe}_complexes/${probe}.confidence_summary.tsv
mv candidate_list.txt ${probe}_complexes
cp ${probe}.combined_plddt.tsv ${probe}_complexes/${probe}.combined_plddt.tsv


tar -czf ${probe}_complexes.tar.gz ${probe}_complexes

##do after moving back to cluster
#mv ${probe}_complexes.tar.gz /staging/s/saodonnell/af3_proteomes/${dataset}
#rm -r ${probe}_complexes



