# Project Description

This is a pipeline I developed for my research which downloads, extracts data from (in parallel), and runs a tanimoto similarity search on every compound available on PubChem.

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

For each of the commands in 'Usage', add the --test flag to the end of the command to run a smaller version of the workflow for testing purposes.

Example:

```bash
python -u download_pubchem_compounds.py --test
```

This will only download the first 3 compound files from pubchem.

## Usage

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

