"""
Extract properties of interest from .sdf.gz files obtained from Pubchem in a directory and save them as individual .zst files.

This workflow supports multiprocessing and can be run in test mode.

for more info, run:
    python extract_data_from_pubchem_sdf.py --help
    
"""

import os
import gzip
import argparse
import multiprocessing as mp
import polars as pl
from tqdm import tqdm
from neattime import neattime
from rdkit import Chem

# properties to extract from pubchem's .sdf files. 
PROPERTIES_TO_EXTRACT_FROM_MOLS = ['PUBCHEM_COMPOUND_CID', 'PUBCHEM_IUPAC_NAME', 'PUBCHEM_SMILES']

NUM_FILES_TO_TEST = 5
NUM_MOLS_TO_TEST = 100


def get_gzipped_sdf_filenames(directory: str) -> list:
    """ Get all filenames in a directory that end with .sdf.gz

    Args:
        directory (str): name of directory to look for .sdf.gz files

    Returns:
        list: list of filenames that end with .sdf.gz
    """
    return [filename for filename in os.listdir(directory) if filename.endswith('sdf.gz')]


def process_gzipped_sdf_file(gzipped_sdf_filepath: str, test: bool = False) -> dict:
    f""" Process a gzipped .sdf file and extract the properties of interest. 

    Properties of interest are defined in PROPERTIES_TO_EXTRACT_FROM_MOLS

    Args:
        gzipped_sdf_filepath (str): path to the gzipped .sdf file
        test (bool, optional): if True, only process the first {NUM_MOLS_TO_TEST} molecules. Defaults to False.

    Returns:
        dict: dictionary of properties of interest and their values (as lists)
    """

    data = {property: [] for property in PROPERTIES_TO_EXTRACT_FROM_MOLS}

    print(f'Processing {gzipped_sdf_filepath}')

    with gzip.open(gzipped_sdf_filepath, 'rb') as f:
        suppl = Chem.ForwardSDMolSupplier(f)

        for i, mol in enumerate(suppl):
            if mol is None:
                continue

            for property in PROPERTIES_TO_EXTRACT_FROM_MOLS:
                try:
                    data[property].append(mol.GetProp(property))
                except KeyError:
                    print(f'Molecule {i} in {gzipped_sdf_filepath}: Property {property} not found in molecule')
                    data[property].append(None)
            
            if test and i == 100:
                break
        
    return data


def save_data_to_zst(data: dict, output_dir: str, filename: str) -> None:
    """Save data in the form of a dictionary to a compressed zst file

    Note: output_dir and filename are combined to form the full path to the file

    Args:
        data (dict): dictionary of data to save
        output_dir (str): directory to save the data
        filename (str): name of the file to save the data as
    """

    df = pl.from_dict(data)

    filepath = os.path.join(output_dir, filename)

    df.write_parquet(file=filepath, compression = 'zstd', compression_level=8)

    print(f'Saved data to {os.path.join(output_dir, filename)}')


def process_gzipped_sdf_file_and_save(gzipped_sdf_filepath: str, output_dir: str, test: bool = False):
    f"""Process a gzipped .sdf file and save the data to a zst file. This is a convenience function needed for multiprocessing.
    It should mimic the behavior of the single process version by combining the process_gzipped_sdf_file and save_data_to_zst functions.

    Args:
        gzipped_sdf_filepath (str): path to the gzipped .sdf file
        output_dir (str): directory to save the data
        test (bool, optional): if True, only process the first {NUM_MOLS_TO_TEST} molecules. Defaults to False.
    """
    data = process_gzipped_sdf_file(gzipped_sdf_filepath, test=test)

    filename = os.path.basename(gzipped_sdf_filepath).replace('.sdf.gz', '.zst')

    save_data_to_zst(data, output_dir, filename)

    print(f'Processed {gzipped_sdf_filepath}')


def cleanup():
    """
    Cleanup function to cancel any paused processes and re-run for any files that failed

    Had a problem where the pipeline got stuck on 2 files and I had to manually kill the process. Not sure what was the problem
    """
    raise NotImplementedError('This function is not implemented yet')


def main(args):
    if args.output_dir is None:
        args.output_dir = f'extracted_pubchem_data/{'TEST_' if args.test else 'full_processed_data_'}{neattime()}/'

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    gzipped_sdf_filenames = get_gzipped_sdf_filenames(args.input_dir)

    # single process

    if args.num_processes == 1:
        for i, gzipped_sdf_filename in enumerate(tqdm(gzipped_sdf_filenames, desc='Processing .sdf.gz files')):
            data = process_gzipped_sdf_file(gzipped_sdf_filepath=os.path.join(args.input_dir, gzipped_sdf_filename), test=args.test)

            out_filename = gzipped_sdf_filename.replace('.sdf.gz', '.zst')

            save_data_to_zst(data, args.output_dir, out_filename)

            if args.test and i == 4:
                break
    
    else:
        if args.num_processes == -1:
            num_processes = mp.cpu_count()
        else:
            num_processes = args.num_processes

        print(f'Using {num_processes} processes')

        # somehow this fixes a problem with the pool getting stuck. See https://pythonspeed.com/articles/python-multiprocessing/
        with mp.Pool(processes=num_processes) as pool:

            gzipped_sdf_filepaths = [os.path.join(args.input_dir, filename) for filename in gzipped_sdf_filenames]

            num_files = len(gzipped_sdf_filepaths)

            pool.starmap(
                process_gzipped_sdf_file_and_save,
                list(zip(gzipped_sdf_filepaths,
                [args.output_dir] * num_files,
                [args.test] * num_files))[:NUM_FILES_TO_TEST if args.test else None]
            )


def parse_args():
    parser = argparse.ArgumentParser(description=f'Extracts properties specified by PROPERTIES_TO_EXTRACT_FROM_MOLS ({PROPERTIES_TO_EXTRACT_FROM_MOLS}) from .sdf.gz files for each molecule obtained from Pubchem in a directory and saves them as individual .zst files')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing pubchem data files as .sdf.gz')
    parser.add_argument('--output_dir', type=str, required=False, default=None, help='Directory to save the extracted data. If not specified, will be saved to extracted_pubchem_data/ with a timestamp')
    parser.add_argument('--num_processes', type=int, default=-1, help='Number of processes to use. -1 means use all available cores')
    parser.add_argument('--test', action='store_true', help=f'Test mode. Only process {NUM_FILES_TO_TEST} files and {NUM_MOLS_TO_TEST} molecules. This is useful for debugging and testing the pipeline. The output will be saved to a directory with TEST_ prefix')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)
    