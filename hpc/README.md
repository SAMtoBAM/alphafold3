# Credits
This pipeline was developped by the **Computational Biology Research Group** at the **GLBRC** (https://www.glbrc.org/). This repository was forked from the original repository (🔗 [link](https://github.com/sameerd/alphafold3/tree/hpc/hpc))  to add clarifications to the instructions. 

# About AlphaFold3
to do

# Cyberinfrastructure Choice
## Is Running AlphaFold3 on CHTC the right option for me?


## Scenarios where running AlphaFold3 on CHTC is a good idea:
- Uncommon ligands
- Large number of novel sequences (e.g. 1000s of jobs)

# Instructions for running AF3 on CHTC

## Pipeline description

This directory has two pipelines for running Alphafold3 on CHTC. It uses an
apptainer container (`.sif` file) (5GB) in two pipelines.

1. `data_pipeline.sh`: The data pipeline searches the genetic databases for matches of interest. It
   copies and decompresses the databases (252GB --> 650GB) to an execute node
   and runs the searches in parallel there. No GPU is needed for this pipeline
2. `inference_pipeline.sh`: The inference pipeline takes the output of the data pipeline and creates a
   structure using a GPU


## **⚠️ Important Notes**

* **Obtain the Alphafold3 model parameters first**: The inference pipeline needs
  [your own copy of the Model
  Parameters](../README.md#obtaining-model-parameters) to run. Google Deepmind
  has terms of use and you need to agree to these, download and save the
  parameters to your own staging directory. Make sure that this file is
  protected by your own file permissions `chmod 600 af3.bin*` and then
  edit `MODEL_PARAM_FILE` variable in the in the [inference pipeline submit script](./inference_pipeline.sub)
* **Batch jobs to avoid using excess CHTC resouces**: Contrary to the usual CHTC
  approach of slicing jobs finely we need to batch jobs as much as possible to
  avoid transferring too much data to/from the staging directory. It is
  possible to put multiple json config files in each of the the input
  directories [job1/data_inputs](./job1/data_inputs/) and
  [job1/inference_inputs](./job1/inference_inputs/) and also create multiple
  `job` directories and run them all at the same time.
* **Avoid copying data**: The apptainer container and the databases in
  `/staging/groups/glbrc_alphafold/af3/` are created and downloaded from public
  sources and are readable by anyone who can access the staging directory. This
  is done so that there is no need to make any copies of this data on the
  staging directory. According to the TERMS of USE, only the model parameters
  need your own copies as you cannot run this software without receiving the
  model parameters directly from Google Deepmind.


## Step-by-step instructions

1. Obtain the model parameters, put it in your own staging directory and
   protect it with your file permissions. `chmod 600 af3.bin.zst`.

> [!WARNING]
> Do not use SQUID since the data should be private, not public.

3. Download the following files to your chtc submit server home directory
   * [data_pipeline.sh](./data_pipeline.sh)
   * [data_pipeline.sub](./data_pipeline.sub)
   * [inference_pipeline.sh](./inference_pipeline.sh)
   * [inference_pipeline.sub](./inference_pipeline.sub)
   * [Makefile](./Makefile) # optional to help run the scripts
4. Create the input directories as a sub directory of the directory that you
   downloaded the above files into.
   ```shell
    mkdir -p job1/data_inputs job1/inference_inputs
   ```
5. Edit the `MODEL_PARAM_FILE` variable in `inference_pipeline.sub` to point
   to the model parameters you got from Google Deepmind in Step 1.

6. Put the config `json` files in the in the [job1/data_inputs/](./job1/data_inputs/)
   directory. An example config to test is available as
   [`fold_input.json`](./test/input/fold_input.json) in
   the [(README.md file for Alphafold3](../README.md)

7. Run a test with the small databases first
   ```shell
   condor_submit USE_SMALL_DB=1 data_pipeline.sub
   ```

   Check and then delete test output files. You will find a file named `2pv7.data_pipeline.tar.gz` located under `job1/2pv7.data_pipeline.tar.gz`

   Unzip the file and take a look at its contents:
   ```
   cd job1
   tar xvf 2pv7.data_pipeline.tar.gz
   cat job1/2pv7_data.json
   ```

You should see this:

```
{
  "dialect": "alphafold3",
  "version": 2,
  "name": "2PV7",
  "sequences": [
    {
      "protein": {
        "id": [
          "A",
          "B"
        ],
        "sequence": "GMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG",
        "modifications": [],
        "unpairedMsa": ">query\nGMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG\n",
        "pairedMsa": ">query\nGMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG\n",
        "templates": []
      }
    }
  ],
  "modelSeeds": [
    1
  ],
  "bondedAtomPairs": null,
  "userCCD": null
}
```

7. If the file `job1/2pv7_data.json` looks ok, proceed to run the `data_pipeline.sub` on the whole dataset.

   ```
   condor_submit data_pipeline.sub
   ```
   
8. This will create a `job1/2pv7.data_pipeline.tar.gz` file.

Note the time stamp, and likely larger file size. Move it to the `job1/inference_inputs` folder:

   ```
   mv job1/*.data_pipeline.tar.gz job1/inference_inputs/.
   ```

9. The run inference pipeline:

   ```
   condor_submit inference_pipeline.sub
   ```

The output will be a file named `.inference_pipeline.tar.gz`. Once unzipped, it should contain the follow files:

```
2pv7.inference_pipeline.tar.gz

# Unzipped
2pv7_confidences.json
2pv7_data.json
2pv7_model.cif
2pv7_summary_confidences.json
ranking_scores.csv
seed-1_sample-0
seed-1_sample-1
seed-1_sample-2
seed-1_sample-3
seed-1_sample-4
TERMS_OF_USE.md
```


## What to do next?

### Run this on multiple sequences

Given a multifasta file, where each is a sequence to be predicted, the script [set_up_directory.py](./set_up_directory.py) can be used to automatically create the folder structure for the data_pipeline.sub to run.
The python scripts takes in 2 arguments: the fasta file, and a batch_size. The batch sized number is used to determine how many `.json` files are created under each `jobN` directory.

1. Download a copy of `set_up_directory.py` and make it executable:

```
wget [...]
chmod +x set_up_directory.py
```

2. Use the python script directly from the access point:
```
python set_up_directory.py test.fasta 10
```

3. Run the data_pipeline.sub script:
```
condor_submit data_pipeline.sub
```

Since the queue statement stays queue directory from job*, it automatically knows to submit a job for each folder started with the naming pattern "job".

4. Move the tar.gz file into the inference_input folders.
A (or multiple) tar.gz file is created under `jobN`, but they all need to be moved to `jobN/inference_inputs/.`

```
for dir in job*/; do mv "$dir"*.tar.gz "$dir"inference_inputs/; done
```

List the directory contents to make sure files have been moved properly:
```
ls job*/inference_inputs/
```

5. Run the inference script for all the samples:

```
condor_submit inference_pipeline.sub
```


### Visualize the model results (`.cif`) files

To visualize the model results (`.cif`) file you will need a program like PyMol.
UW-Madison researchers can go to software.wisc.edu, login with their netID and password, and search for PyMol. In the download options, you can download the software (for your Windows, Mac or Linux) computer *and* there is a link to the license ("PyMOL License File - 3.0 Version or Later"). 

Download and Install Pymol.
Choose your license (the `.lic`) file you just downloaded.

Get a copy of the `cif` file off of CHTC:

In a new terminal window type
```
cd ~/Downloads
sftp netid@ap2002.chtc.wisc.edu
# navigate to the alphafold3 output folder
# For example:
cd ~/AF3/job1
get 2pv7_model.cif
exit
```
The cif file will be downloaded to your Downloads folder on your computer.

Open PyMol. Click on File > Open > Choose the cif file.

## Notes on recreating steps used to create these instructions

These notes below are only needed in case someone wants to recreate building
the apptainer container or the databases


### Creating the apptainer container
The latest version of the apptainer container is stored at
`/staging/groups/glbrc_alphafold/af3/alphafold3.sif` and the commands below are
just in case anyone needs to build it for themselves. This container is built
using the Google Deepmind Alphafold3 publically available code and it's use
should be governed by the same terms of use.

```shell
# Build the container on your own machine (if you want to)
# Docker needs around 11GB of memory to run
# On an 8GB machine, enable swapspace of 4GB so that the `build_data` command
# does not get killed
sudo docker build -t alphafold3 -f docker/Dockerfile .
# apptainer needs a lot less memory but needs around 30GB of disk space
# mostly in it's tmp directory. Can be set with APPTAINER_CACHE and TMPDIR
sudo apptainer build alphafold3.sif docker-daemon://alphafold3:latest


# creating smaller singularity image without build components
cd hpc
sudo docker build -t alphafold3_minimal -f Dockerfile .
sudo apptainer build alphafold3.minimal.sif \
        docker-daemon://alphafold3_minimal:latest
```



### Creating small databases for testing

```
# conda env utils has zstd so activate it
eval "$(micromamba shell hook --shell bash)"
micromamba activate utils


compressed_db_dir= /staging/groups/glbrc_alphafold/af3/db

small_db_dir=/staging/groups/glbrc_alphafold/af3/db_small
mkdir -p ${small_db_dir}
touch "${small_db_dir}/test.zst"
/bin/rm "${small_db_dir}"/*.zst

pushd "$compressed_db_dir"
for archive in *.zst
do
  echo processing "${archive}"
  zstdcat -d "${archive}" | ~/bin/subsample-fasta -n 1000 | zstd - -o "${small_db_dir}/${archive}"
done
popd



```
