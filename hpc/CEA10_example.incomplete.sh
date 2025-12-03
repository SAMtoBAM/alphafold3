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

#############################################
######### 1a. GET AF3 MODEL WEIGHTS #########
#############################################


##the weights for the alphafold3 model are not open source therefore:
##go to this website to ask for permission for the alphafold3 models https://github.com/SAMtoBAM/alphafold3/blob/hpc/README.md#obtaining-model-parameters

## once you have a link to the file, download it, place it is your staging as so /staging/${user}/af3/weights/af3.bin.zst
##then modify the file permissions
chmod 600 /staging/${user}/af3/weights/af3.bin.zst


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
mv af3complex.sif /staging/${user}/af3complex.sif


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

##modify the paths to your singularity image (still assuming it is here: '/staging/${user}/af3complex.sif')
sed -i "s/USERNAME/${user}/" *.complex.sub

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
proteome="../proteins/SAMN28487501.proteins.t1.fa"

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

##if some jobs are not still running AND have no *.tar.gz files these may have not completed for some weird reason
##You can re-run any specific jobs by changing the last line of the submission script 'queue directory matching job*'
##change job* to the job with the missing output e.g. 'queue directory matching job99' and do this for each job without output


##move all msas into a single folder and clean up jobs
mkdir proteome_msa
mv job*/*.data_pipeline.tar.gz proteome_msa/
rm -r jobs/


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
probe="g1367"
##variable for the protein searched against others (can use to simplify summary output files)
probename=LaeA

##have new python script that should combine the probe with all partners in the combined proteom-msa folder
##this will recreate the jobN folder and inference_input folder and placed the combo msa jsons there
python msa_pairing.py \
  --msa_dir ./proteome_msa/ \
  --probe ./proteome_msa/${probe}_data.json \
  --out_dir ./ \
  --batch_size 100


##and now we submit jobs for the inference step
condor_submit inference_pipeline.complex.sub

##if everything completed you should once again see some files in the job*/ folders called '*.inference_pipeline.tar.gz'
##these folders now contain all the data for the structures, for each model and for each sample for each model (currently 5 samples for each of the 20 models)

#############################################
######### 5. STRUCTURE STATISTICS ##########
#############################################

##Assuming everything ran well we now want to evaluate the best model-sample combination per complex

