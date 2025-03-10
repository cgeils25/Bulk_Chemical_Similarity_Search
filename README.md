# Project Description

This is a data pipeline I developed for my research which downloads, extracts data from, and runs a tanimoto similarity search on every compound available on PubChem.

I tried my best to optimize this with parallelization and file compression to reduce storage requirements, but you could probably push things further by writing in C/C++ and with some more clever tricks. I skipped this because I had to publish my paper eventually...

## Setup

The environment for this project is managed with conda. If you don't have conda installed, [view the instructions to download miniconda first](https://www.anaconda.com/docs/getting-started/miniconda/install).

Once conda is installed, run the following to build a suitable environment:

```bash
conda env create --file requirements.yml
```

Next, activate the environment with:

```bash
conda activate bss
```

## Tests

For each of the commands in 'Usage' (aside from `source run_full_pipeline.sh`), add the --test flag to the end of the command to run a smaller version of the workflow for testing purposes.

Example:

```bash
python -u download_pubchem_compounds.py --test
```

This will only download the first 3 compound files from pubchem.

## Usage

### To Run the Entire Pipeline at Once

Run: 

```bash
source run_full_pipeline.sh
```

TODO add some args to this script so users can specify their dataset, threshold, number of cores, etc.

### To Run it in Stages

To download every compound available through pubchem as an SDF file, run:

```bash
python -u download_pubchem_compounds.py
```

(the -u flag is only important if you're writing to a log file. Without it, python just prints everything once the process finishes rather than continuously)

FYI this ended up being 111 GB of compressed (.gz) files for me. Uncompressed, I calculated it would be 909.8 GB, although my implementation leaves the files themselves compressed while processing them.

Next, to extract relevant data (SMILES, pubchem Id, etc.) from the resulting sdf files with parallelization, run:

```bash
python -u extract_data_from_pubchem_sdf.py --input_dir ... --num_processes -1
```

Where '...' is the path to the directory containing the output of `download_pubchem_compounds.py`

Next, to compute the pairwise tanimoto similarity between compounds in your dataset and every compound on pubchem, run:

```bash
python -u compute_tanimoto_similarity.py \
    --comparison_dataset test_comparison_dataset.csv \
    --extracted_pubchem_data_dir ... \
    --num_processes 8
```

Where '...' is the path to the directory containing the output of `extract_data_from_pubchem_sdf.py`

You can replace `test_comparison_dataset.csv` with your data. It should contain a column named "smiles" whose rows are SMILES strings. 

NOTE: be careful with the number of processes for `compute_tanimoto_similarity.py`. When running this on an HPC node with a 48-core intel xeon and 196 GB of RAM, the highest I could do while remaining stable was 16. This is the most memory-intense part of the pipeline so adjust to your needs.  

Finally, filter out compounds without at least one tanimoto similarity score above a given threshold (here 0.8), run:

```bash
python filter_tanimoto_results.py \
    --tanimoto_directory ... \
    --threshold 0.8
```

Where '...' is the path to the directory containing the output of `tanimoto_similarity_search.py`

This last step is of course optional and you could really do whatever you want with the computed tanimoto scores
