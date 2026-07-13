####These scripts were used in order to analyse the protein-protein binding of proteins of interest against entire proteomes
##The script does so by employing AF3Complex; a version of alphafold3 that has been modified to improve complex detection
##the idea behind this development was to screen proteins of interest against protein databases in order to see what interactions they have
##by doing this we may be able to learn about the protein function, most useful within proteins with no known function.

##the basic pipeline is:
#1. Set up AF3Complex
#2. Run AF3Complex multiple sequence alignment on the entire proteome of your genome
#3. Select a single protein of interest as a probe and combine this with the entire proteome
#4. Run the protein structural inference and complex confidence
#5. Extract the confidence scores and isolate candidates from background

#########################################################################
############################## 1. Set up ################################
#########################################################################

##All steps are laid out here: https://github.com/SAMtoBAM/alphafold3/tree/hpc/hpc

##put your username here (will be used a few times in the set up scripts)
user="saodonnell"

##need to set this variable as the first letter of your user
##this is now in the directory systme for users in staging (a subdirectory split by first letters in the user name)
firstletter=$( echo ${user} | cut -c1-1 )

#############################################
######### 1a. GET AF3 MODEL WEIGHTS #########
#############################################


##the weights for the alphafold3 model are not open source therefore:
##go to this website to ask for permission for the alphafold3 models https://github.com/SAMtoBAM/alphafold3/blob/hpc/README.md#obtaining-model-parameters

## once you have a link to the file, download it, place it is your staging as so /staging/${firstletter}/${user}/af3/weights/af3.bin.zst
##then modify the file permissions
chmod 600 /staging/${firstletter}/${user}/af3/weights/af3.bin.zst


#############################################
########### 1b. CONDA ENVIRONMENT ###########
#############################################

##conda environment used while setting up (outside of the singularity container used to actually run alphafold)
##for now just need biopython installed in-order to run the python script that sets up the directory/scripts for each run

##first get miniconda installed in /home/USERNAME/
#cd ~
#mkdir -p ~/miniconda3
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
#bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
#rm ~/miniconda3/miniconda.sh
#source ~/miniconda3/bin/activate
#conda init --all
##exit out of terminal then reconnect to intialise

##now create actual environment
#conda create -n AF3 conda-forge::biopython

conda activate AF3

#############################################
######### 1c. GET SINGULARITY IMAGE #########
#############################################


##a docker image for AF3Complex has already been created so you just need to pull the docker image and create a singularity one
##using singularity/apptainer to avoid permission issues with docker
singularity pull af3complex.sif docker://samtobam/af3complex:af3
##move this singularity image to staging (this is where we will expect the image to be during submissions)
mv af3complex.sif /staging/${firstletter}/${user}/af3complex.sif


#############################################
########### 1d. PROTEIN DATABASE ############
#############################################

##This has already been downloaded and stored in an accessible region of the CHTC
##this takes up ~250GB compressed and ~600GB decompressed therefore it is nice to have stored somewhere
#/staging/groups/glbrc_alphafold/af3

#########################################################################
########################## 2. AF3complex MSAs ###########################
#########################################################################

##########
########## FROM NOW ALL ANALYSES ARE TO BE TAILORED TO DIFFERENT GENOMES; THEREFORE RERUNNING THE MSA ANALYSIS OF ENTIRE PROTEOMES
##########

##give your dataset a name (used as output folder for all results); can just be the name of the genome/where you proteome comes from
dataset=CEA10_SAMN28487501

cd ~/alphafold3/
mkdir ${dataset}
cd ${dataset}

##folder to store all the confidence summaries for each probe for easy access (without taking up much space)
mkdir confidence_summaries/

#############################################
#### 2a. Download analysis scripts etc ######
#############################################

##All the scripts required for running this analysis have already been written and designed for the CHTC


##submission scripts and their associated bash scripts for running AF3Complex on the cluster have already been created
##just need to download them
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/inference_pipeline.complex.sh
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/inference_pipeline.complex.sub
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/data_pipeline.complex.sub
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/data_pipeline.complex.sh

##modify the paths to your singularity image (still assuming it is here: '/staging/${firstletter}/${user}/af3complex.sif')
sed -i "s/USERNAME/${user}/" *.complex.sub
sed -i "s/FIRSTLETTER/${firstletter}/" *.complex.sub

##we also need a python script that will organise the input proteins
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/set_up_directory.py
##change the permissions
chmod +x set_up_directory.py

##a script that will combine the MSAs using a defined probe
wget https://raw.githubusercontent.com/SAMtoBAM/alphafold3/hpc/hpc/AF3complex/msa_pairing.py
chmod +x msa_pairing.py




#########################################################################
############################ CEA10 proteome #############################
#########################################################################


#############################################
########## 2b. Create input files ###########
#############################################

##for this example I have placed the proteome fasta files in a folder called proteins
## SAMN28487501.proteins.t1.fa (proteome of CEA10; annotated in-house and extracting only the protein associated with the first transcript of each gene)

##now we just need a fasta file with all the proteins of interest 
##each one will get an individual MSA which can be used to predict stucture and complexes over and over again
##THIS DATASET CONTAIN SPLIT PROTEINS BASED ON PROTEIN DOMAINS DUE TO THEM BEING TOO LARGE FOR INFERENCE (see below)
proteome="../proteins/SAMN28487501.proteins.t1.split.fa"

##first get list of proteins from proteome file
##will be used when generating the summary statistics file
grep '>' "${proteome}" | sed 's/>//g' > full_protein_list.txt

##just need one more variable; the number of complexes to combine into single jobs that will be submitted
##for now 30 seems to work well (balancing the time taken by the first step to decompress the protein database and how many proteins will go over the initial memory request...its hard to compute)
##this will create a jobX with each job having subdirectories, one 'af_input' containing a json file with information for ${proteinsperjob} number of complexes
proteinsperjob="30"

##make sure the conda env is activate 'conda activate AF3'
##run the set up script which generates single protein jsons accepted by alphafold3 for every protein in your proteome
python set_up_directory.py ${proteome} ${proteinsperjob}


#############################################
###### 2c. MULTIPLE SEQUENCE ALIGNMENT ######
#############################################

##now we have our input json per job
##we can launch a batch submission for all jobs that performs the first step in the analysis, multiple sequence alignment
##this queries each protein against the protein database already on the cluster and aligns them then spits out a json, per complex, with the alignment information
##this step takes the longest and the most chances for going over the requested memory

condor_submit data_pipeline.complex.sub

##once the alignment step is complete ou should see in job*/ that there are a number of *.data_pipeline.tar.gz folders
##each one will contain a json file with the complex information and the alignments
##make sure all are there before continuining, there should be the same number of *.data_pipeline.tar.gz files as proteins in your proteome

##can check by running this
outputcount=$( ls job*/*.data_pipeline.tar.gz | wc -l )
inputcount=$( grep ^'>' ${proteome} | wc -l )

if [[ "$outputcount" -eq "$inputcount" ]]; then
  echo "Looks good! Number of output files match input (${inputcount})"
else
  echo "Number of proteome input sequences does not match number of output sequences"
  echo "Output folders (XXXX.data_pipeline.tar.gz) count: $outputcount"
  echo "Proteome protein count: $inputcount"
  echo "Are some jobs still running? Did some jobs not complete?"
	echo "Checking for job directories without .data_pipeline.tar.gz files..."
	for d in job*/
	do
		[[ -d "$d" ]] || continue
		if ! ls "$d"/*.data_pipeline.tar.gz >/dev/null 2>&1
		then
			echo "$d has no output files (XXXX.data_pipeline.tar.gz)"
		fi
	done
fi


##several large proteins g939, g2318, g2671, g3537, g5851, g6902, g7549, g8094, g8435 and g9331 were split into two
##in both cases the proteins were split in between functional domains/structural features
##each are now two seperate proteins sequences -1 and -2 e.g. g9331-1 and g9331-2
##for g8094 and g9331, they were split because they too much RAM during MSA
##the reamining all made it through the MSA process intitially but require too much GPU memory (VRAM) during inference


##if some jobs are not still running AND have no *.tar.gz files these may have not completed for some weird reason
##You can re-run any specific jobs by changing the last line of the submission script 'queue directory matching job*'
##change job* to the job with the missing output e.g. 'queue directory matching job99' and do this for each job without output


##move all msas into a single folder and clean up jobs
mkdir proteome_msa
mv job*/*.data_pipeline.tar.gz proteome_msa/
rm -r job*/

##just zip up the proteome_msa folder as it takes up a lot of space for the mean time
##sadly will have to unzip it each time wanting to run another probe (but that is quite fast to do)
tar -czf proteome_msa.tar.gz proteome_msa/
mkdir /staging/${firstletter}/${user}/af3_proteomes/
mkdir /staging/${firstletter}/${user}/af3_proteomes/${dataset}
mv proteome_msa.tar.gz /staging/${firstletter}/${user}/af3_proteomes/${dataset}/
rm -r proteome_msa/


########################################################
########## 3. AF3Complex STRUCTURAL INFERENCE ##########
########################################################

##########
########## FROM NOW ALL ANALYSES ARE TO BE TAILORED TO DIFFERENT PROBES; THEREFORE RERUNNING THE INFERENCE ANALYSIS OF A SINGLE PROBE VS THE WHOLE PROTEOME
##########


##this section will now combine each of the MSAs with a 'probe' protein
##therefore inference will be looking at the interaction of the predicted protein structures for all proteins in the proteome with your probe protein


##Deviations from defaults runs of AF3Complex
	## 1. The inference of each structure tests 20 different models
		##This is to find even better models for each complex and aligns with the paper showing that more models helps find better structures even up until 1000 but where 20 shows huge improvements then decreases
	## 2. Assumes you are interested in pairwise protein-protein interactions
	## 3. Uses the alphafoldserver json input format for MSA


##now we just need to know the probe
##the 'probe' file should contain a single protein which will be combined against all proteins in the proteome as complexes
##the number of complexes analysed will therefore be the number of proteins in the 'proteome'
##this name should be the name of the protein in the proteome file and therefore have an MSA file called '${probe}_data.json' 
probe="g3045"

##if already compressed (just send the decompressed files to here)
tar -xzf  /staging/${firstletter}/${user}/af3_proteomes/${dataset}/proteome_msa.tar.gz -C ./

##have new python script that should combine the probe with all partners in the combined proteom-msa folder
##this will recreate the jobN folder and inference_input folder and placed the combo msa jsons there
python msa_pairing.py \
  --msa_dir ./proteome_msa/ \
  --probe ./proteome_msa/${probe}.data_pipeline.tar.gz \
  --out_dir ./ \
  --batch_size 10

##this step can take awhile as it needs to open up and then join the protein sequences
##if it doesn't finish generating all files (before timing out for example) just run it again and it'll skip those already completed (you can also remove the last job folder created just to be safe)


##and now we submit jobs for the inference step
condor_submit inference_pipeline.complex.sub

##if everything completed you should once again see some files in the job*/ folders called '*.inference_pipeline.tar.gz'
##these folders now contain all the data for the structures, for each model and for each sample for each model (currently 5 samples for each of the 20 models)


####### FIRST CHECK THAT ALL BATCH JOBS FINISHED

##automated check that everything worked
outputcount=$( ls job*/*.inference_pipeline.tar.gz | wc -l )
inputcount=$( ls job*/inference_inputs/*.data_pipeline.tar.gz | wc -l )

if [[ "$outputcount" -eq "$inputcount" ]]; then
  echo "Looks good! Number of output files match input (${inputcount})"
  echo "Next check if all folders actually contain model outputs (sometimes an error for a single complex occurs, the run doesn't finish but the output folder is still generated)"
else
  echo "Number of input inference files does not match number of output folders"
  echo "Output folders (job*/*.inference_pipeline.tar.gz) count: $outputcount"
  echo "Input folders count (job*/inference_inputs/*.data_pipeline.tar.gz): $inputcount"
  echo "Are some jobs still running? Did some jobs not complete?"
fi

##if all looks good and everything ran then move all the outputs to another folder
mkdir proteome_inference
mv job*/*inference_pipeline.tar.gz proteome_inference/


####### SECOND CHECK THAT ALL COMPLEXES WERE INFERRED

##check if any of the interence output folders are actually missing the output
##generally jobs may not finish because they went over the VRAM limit of the GPU used
##but not to worry, we can resubmit requesting a higher VRAM GPU

##just look for the model *.cif file in the output folder
##check if empty and if so; copy the input into another folder for rerunning the inference
mkdir rerun-job
mkdir rerun-job/inference_inputs
missing_cif="0"
for f in proteome_inference/*.inference_pipeline.tar.gz; do
    [[ -e "$f" ]] || continue

    if ! tar -tf "$f" | grep -qi '\.cif$'; then
        echo "WARNING: no .cif found in $f"
        ((missing_cif++))

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
echo "The new file has a higher VRAM minimum (from 16GB to 60GB) and only runs on only the data in the folder 'job-rerun'"
cat inference_pipeline.complex.sub | sed 's/CUDAGlobalMemoryMb > 48000/CUDAGlobalMemoryMb > 80000/' | sed 's/job\*/rerun-job/' > inference_pipeline.complex.rerun.sub
fi

##resubmit the failed inferences with the higher VRAM minimum request
condor_submit inference_pipeline.complex.rerun.sub
##now move the output into the proteome_inference folder again (overwriting the old empty ones)
mv rerun-job/*inference_pipeline.tar.gz proteome_inference/

##now run the check above again!
##if this didn't work, likely you will need to split the protein (trying not to break domains); re-run the MSA step; then rerun inference. or just skip it.


##once finished, tidy up the input data
rm -r job*
rm -r rerun-job


#############################################
########## 4. STRUCTURE STATISTICS ##########
#############################################


##IF ONLY ONE SEED WAS GIVEN (DEFAULT) AND GOOD FOR A FIRST BIG PROTEOME WIDE SEARCH TO JUST DO A SINGLE SEED

##can now clear out the inference input folder
##we can remove all the different analyses (keeping just the best one as per the last seed run) as each one can contain a lot of unnecceessary data
##will also remove the MSA used for inference since that is already stored in the individual MSAs (the largest file in the output) (also makes it much quicker to open/compress etc)
##need to be a bit prudent as to what is kept here as the data will really begin to pile up with the all-v-all for 10k proteins

##ONLY IF THERE IS NOT SEVERAL SEEDS THAT WERE RUN!

mkdir proteome_inference_clean
for f in proteome_inference/*.inference_pipeline.tar.gz
do
  complex=$(basename "${f}" .inference_pipeline.tar.gz)
  if [[ ! -f "proteome_inference_clean/${complex}.inference_pipeline.tar.gz" ]]
  then
  echo "Processing $complex"
  tmpdir=$(mktemp -d)
  # Extract only kept files into temp dir
  tar -xzf "${f}" --exclude='seed*' --exclude='*_data.json' -C "$tmpdir"
  # Repack
  tar -czf "proteome_inference_clean/${complex}.inference_pipeline.tar.gz" -C "$tmpdir" .
  rm -rf "$tmpdir"
  fi
done

rm -r proteome_inference
mv proteome_inference_clean proteome_inference

##Assuming everything ran well we now want to evaluate the best model-sample combination per comple
##already the best seed has been evaluated for the last model/seed run. So if only one seed was given (default) then run the below getting information on the best seed for that seed
##this handles getting confidence scores whether or not the file is completed and fills in NAs if nothing
##folder search also needs to allow for some letter being made lowercase by alphafold

##use list of proteins to analyse
echo "protein1;protein2;fraction_disordered;has_clash;pLDDT;ptm;iptm;ranking_score" | tr ';' '\t' > confidence_summary.tsv
cat full_protein_list.txt | while read -r pair
do
    folder=$(find proteome_inference \
    -maxdepth 1 \
    -type f \
    -iname "${probe}_${pair}.inference_pipeline.tar.gz" \
    -print -quit)
    protein=$(basename "$folder" .inference_pipeline.tar.gz)

    ## defaults
    pLDDT="NA"
    iptm="NA"
    ptm="NA"
    RS="NA"
    FD="NA"
    HC="NA"

    ## pLDDT
    cif_out=$(tar -axf "$folder" "./${protein}_model.cif" -O 2>/dev/null)
    if [[ -n "$cif_out" ]]; then
        pLDDT=$(echo "$cif_out" \
            | awk '/_ma_qa_metric_global.metric_value/ {print $2; exit}')
    fi

    ## JSON metrics (SAFE)
    if tar -tf "$folder" "./${protein}_summary_confidences.json" &>/dev/null; then
        read iptm ptm RS FD HC < <(
            tar -axf "$folder" "./${protein}_summary_confidences.json" -O \
            | jq -r '
                [
                  .iptm,
                  .ptm,
                  .ranking_score,
                  .fraction_disordered,
                  .has_clash
                ] | map(. // "NA") | @tsv
            '
        )
    fi

    ## output
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$probe" "$pair" "$FD" "$HC" "$pLDDT" "$ptm" "$iptm" "$RS"

done >> confidence_summary.tsv


##the below is another option that can handle if more seeds are provided/used for inference 
##therefore extracting all the relevant scores (ptm, iptm, ranking score)
##but first we need to pick, for each model, the best sample based on the ranking score and just keep all the stats for that one
##and for the first summary output we can output all the stats per model

#echo "complex;seed;best_sample;fraction_disordered;has_clash;pLDDT;ptm;iptm;ranking_score" | tr ';' '\t' > confidence_summary.all_seeds.tsv
#ls job*/*.tar.gz | while read folder
#do
###individually run through each seed and get the best results per sample (based on highest confidence score)
#tar -tf $folder | awk -F "_" '{print $1}' | sort -u | grep seed | awk -F "/" '{print $2}' | while read seed
#do
#tar -tf $folder | grep "${seed}_" | awk -F "_" '{print $2}' | awk -F "/" '{print $1}' | sort -u | while read sample
#do
#tar --wildcards  -axf $folder ./${seed}_${sample}/summary_confidences.json -O | sed 's|"||g' | awk -F ":" -v protein="$protein" -v seed="$seed" -v sample="$sample" '{if($1 == " ranking_score") {print seed"\t"sample"\t"$2}}' | sed 's/,//g'
#done | sort -k3n | tail -n1 | awk '{print $2}' | while read bestsample
#do
###get name of protein used
#protein=$( echo $folder | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
###read the compressed output folder and extract the pLDDT scores from the *_model.cif file
#pLDDT=$( tar -axf $folder ./${seed}_${bestsample}/model.cif -O | grep _ma_qa_metric_global.metric_value | awk -F " " '{print $2}' )
###read the compressed output folder and extract the iptm, ptm and ranking_scores from the confidence summary file
###then print all the important stats for each protein/complex
#tar --wildcards  -axf $folder ./${seed}_${bestsample}/summary_confidences.json -O | sed 's|"||g' | awk -F ":" -v protein="$protein" -v pLDDT="$pLDDT" -v seed="$seed" -v bestsample="$bestsample" '{if($1 == " iptm") {iptm = $2} else if($1 == " ptm") {ptm = $2} else if($1 == " ranking_score") {RS = $2} else if($1 == " fraction_disordered") {FD = $2} else if($1 == " has_clash") {HC = $2}} END{print protein"\t"seed"\t"bestsample"\t"FD"\t"HC"\t"pLDDT"\t"ptm"\t"iptm"\t"RS}' | sed 's/,//g'
#done
#done 
#done >> confidence_summary.all_seeds.tsv

###now get the best model per sample/protein-complex based on the summary_score again
#awk 'NR==1 {print; next} {
#  key=$1
#  val=$NF
#  if (val > max[key]) {
#    max[key]=val
#    line[key]=$0
#  }
#}
#END {
#  for (k in line) print line[k]
#}' confidence_summary.all_seeds.tsv > confidence_summary.best_seed.tsv



##now we could just evaluate this manually, but what about if we have compared each probe against many example
##perhaps we only care about the strong outliers?
##to do this we can calculate emperical p-values and select everything over 0.005

##use emperical p-values to find the candidates k-mers
##Step 1: Count total number of scores
#total=$( tail -n+2 confidence_summary.tsv  | wc -l)

##Step 2: Compute empirical p-values and label 
#cat confidence_summary.tsv | awk -F "\t" '{print $1"\t"$7}' | sort -k2,2n -t$'\t' | awk '
#NR==1 {print $0 "\tp_empirical\tcandidates"; next}
#{
#    scores[NR-1] = $2
#    lines[NR-1] = $0
#    n = NR-1
#}
#END {
#    for(i=1;i<=n;i++){
#        count=0
#        for(j=1;j<=n;j++){
#            if(scores[j] <= scores[i]) count++
#        }
#        p = 1 - (count/n)
#        label = (p < 0.005 ? "candidates" : "background")
#        print lines[i] "\t" p "\t" label
#    }
#}' | grep candidates | cut -f1 > confidence_summary.candidates_list.txt
##now label the original summary file
#awk 'NR==FNR {c[$1]=1; next} 
#NR==1 {print $0 "\tcandidate_status"; next} 
#{
#  label = ($1 in c ? "candidates" : "background")
#  print $0 "\t" label
#}' confidence_summary.candidates_list.txt confidence_summary.tsv > confidence_summary.candidates.tsv


##now we can use these confidence scores to try and remove some really unlikely complexes in order to save room in storage
##saving all proteins = ~10GB per probe...therefore we would have 10k*10GB=100TB
##so if we remove ~95% of bad calls we can reduce this to around 5TB...more reasonable 

##use emperical p-values to find the good candidates (also use the ranking score and the iptm seperately and take the nonredundant list of candidates)
##with a cut of of 0.5 (lenient here in order to be forgiving but save ALOT of space)

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
        label = (p < 0.05 ? "candidates" : "background")
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
        label = (p < 0.05 ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 4: get a nonredundant list
cat temp_list.txt | sort -u > candidate_list.txt
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

##get rid of the rest and rename
rm -r proteome_inference
mv proteome_inference_candidates proteome_inference

##reduce the file numbers for storage by zipping up all the already zipped results
tar -czf proteome_inference.tar.gz proteome_inference

cp confidence_summary.tsv confidence_summaries/${probe}.confidence_summary.tsv

##pack up the output for the probe and move it off for temp storage
mkdir ${probe}_complexes
mv proteome_inference.tar.gz ${probe}_complexes
mv confidence_summary.tsv ${probe}_complexes
mv candidate_list.txt ${probe}_complexes

tar -czf ${probe}_complexes.tar.gz ${probe}_complexes
mv ${probe}_complexes.tar.gz /staging/${firstletter}/${user}/af3_proteomes/${dataset}
rm -r proteome_inference
rm -r ${probe}_complexes





#########################################################################
############################ 5. XXXXXXXXXX ##############################
#########################################################################


###temp scripts for analysing some results in terms of overlapping good candidates
cd confidence_summaries

##use the meperical p-values to get good candidates (however now use a higher cut-off e.g. 0.01)
cutoff=0.01

ls *.confidence_summary.tsv | while read file
do
probe=$( echo $file | sed 's/.confidence_summary.tsv//g' )


##Step 1: Count total number of scores
total=$( tail -n+2 ${file} | wc -l)

##Step 2: Compute empirical p-values using the ranking score and grab those below a cut off pval (and ignoring NAs)
tail -n+2 "${file}" |
awk -F"\t" -v cutoff="$cutoff" '
$8 != "NA" {
    scores[++m] = $8
}
{
    lines[NR] = $0
    vals[NR]  = $8
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
        label = (p < cutoff ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2 > temp_list.txt

##Step 3: same as above but using only the iptm
tail -n+2 ${file} |
awk -F"\t" -v cutoff="$cutoff" '
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
        label = (p < cutoff ? "candidates" : "background")
        print lines[i] "\t" p "\t" label
    }
}' | grep -P "\tcandidates$" | cut -f2  >> temp_list.txt

##Step 4: get a nonredundant list
cat temp_list.txt | sort -u > ${probe}.candidate_list.${cutoff}.txt
rm temp_list.txt

done

##get a list of all candidates and which probes they are associated with
echo "candidate;probe" | tr ';' '\t' > candidate_list.${cutoff}.overlaps.tsv
awk '
{
  files[$0] = files[$0] ? files[$0] "," FILENAME : FILENAME
}
END {
  for (s in files)
    print s "\t" files[s]
}
' *.candidate_list.${cutoff}.txt | sed "s/.candidate_list.${cutoff}.txt//g" >> candidate_list.${cutoff}.overlaps.tsv


## get the total number of overlapping candidates for each probe combo
echo "probes;count" | tr ';' '\t' > candidate_list.${cutoff}.overlaps.upset_combinations.tsv
tail -n+2 candidate_list.${cutoff}.overlaps.tsv | awk -F'\t' '
{
  n = split($2, f, ",")
  asort(f)
  combo = f[1]
  for (i=2; i<=n; i++) combo = combo "," f[i]
  count[combo]++
}
END {
  for (c in count)
    print c "\t" count[c]
}
' | sort -k2,2nr >> candidate_list.${cutoff}.overlaps.upset_combinations.tsv

##get a matrix per candidate, presence or absence as candidate per probe
tail -n+2 candidate_list.${cutoff}.overlaps.tsv  | awk -F'\t' '
{
  n = split($2, f, ",")
  for (i=1; i<=n; i++) {
    present[$1][f[i]] = 1
    files[f[i]] = 1
  }
  strings[$1] = 1
}
END {
  # header
  printf "string"
  for (file in files) printf "\t%s", file
  print ""

  # rows
  for (s in strings) {
    printf "%s", s
    for (file in files)
      printf "\t%d", (present[s][file] ? 1 : 0)
    print ""
  }
}
' > candidate_list.${cutoff}.overlaps.upset_matrix.tsv
