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

## Usage

To download every compound available through pubchem as an SDF file, run:

```bash
python -u download_pubchem_compounds.py
```

(the -u flag is only important if you're writing to a log file. Without it, python just prints everything once the process finishes rather than continuously)

FYI this ended up being 111 GB of .gz files for me. 
