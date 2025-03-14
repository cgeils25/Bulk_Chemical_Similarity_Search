# Project Description

This is a data pipeline I developed for my research which downloads, extracts data from, and runs a tanimoto similarity search on every compound available on PubChem.

I tried my best to optimize this with parallelization and file compression to reduce storage requirements, but you could probably push things further by writing in C/C++ and with some more clever tricks. I skipped this because I had to publish my paper eventually...

## Setup

### With Docker

For reproducibility, you can build a docker container for this project. 

### Without Docker

The environment for this project is managed with uv. If you don't have uv installed, [view the instructions for installing uv.](https://docs.astral.sh/uv/getting-started/installation/)

Once uv is installed, build a suitable environment with:

```bash
uv sync
```

## Tests

For each of the commands in 'Usage' (aside from `source run_full_pipeline.sh`), add the --test flag to the end of the command to run a smaller version of the workflow for testing purposes.

Example:

```bash
uv run python -u download_pubchem_compounds.py --test
```

This will only download the first 3 compound files from pubchem.

## Usage

### To Run the Pipeline in Stages

There are 4 components, which I ran one after another. You could string these together in one shell script by specifying input and output directory names, although I wouldn't recommend it because you really only want to do steps 1 and 2 one time. 

1. `download_pubchem_compounds.py` - downloads every compound available on PubChem (119M as of March 2025) as a collection of .sdf.gz files
2. `extract_data_from_pubchem_sdf.py` - extracts relevant information for each molecule in PubChem .sdf.gz files and saves as parquet (.zst) files. See PROPERTIES_TO_EXTRACT_FROM_MOLS to add / remove properties. The most important is PUBCHEM_SMILES which is used to generate morgan fingerprints
3. `compute_tanimoto_similarity.py` - computes pairwise tanimoto similarity between every compound 
4. `filter_tanimoto_results.py` - filter for PubChem compounds with at least one Tanimoto similarity score above a given threshold. 

To view instructions and information about command line arguments for each file, add the '--help' flag. 

Example:

```bash
uv run python download_pubchem_compounds.py --help
```

Output:

```
usage: download_pubchem_compounds.py [-h] [--test] [--output_dir OUTPUT_DIR]

Download all PubChem compounds as .sdf.gz files from the pubchem ftp site (https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/)

options:
  -h, --help            show this help message and exit
  --test                Run in test mode. Only download 3 files.
  --output_dir OUTPUT_DIR
                        Directory to save the downloaded files. If not specified, will be saved to pubchem_data/ with a timestamp.
```
