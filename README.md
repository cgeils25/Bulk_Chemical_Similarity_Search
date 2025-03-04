## Setup

```bash
conda env create --file requirements.yml
```

```bash
conda activate bss
```

## Tests

To download a small sample of pubchem's compounds for testing purposes, run:

```bash
python download_pubchem_compounds.py --test
```

## Usage

To download every compound available through pubchem, run:

```bash
python download_pubchem_compounds.py
```
