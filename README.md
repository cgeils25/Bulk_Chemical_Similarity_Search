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

To download a small sample of pubchem's compounds for testing purposes, run:

```bash
python download_pubchem_compounds.py --test
```

Next, test extracting relevant data (SMILES, pubchem Id, etc.) from the resulting sdf files with:

```bash
python extract_data_from_pubchem_sdf.py --input_dir ... --test
```

Where '...' is the path to the directory containin the output of `download_pubchem_compounds.py`


## Usage

To download every compound available through pubchem as an SDF file, run:

```bash
python -u download_pubchem_compounds.py
```

(the -u flag is only important if you're writing to a log file. Without it, python just prints everything once the process finishes rather than continuously)

FYI this ended up being 111 GB of compressed (.gz) files for me. Uncompressed, I calculated it would be 909.8 GB
