"""
--what does this file do--

To be clear, you can just download the .smi files directly from PubChem, but I wanted to have all the data in one place in case I need it for the future
"""
import os
import gzip
import argparse
import multiprocessing as mp
from tqdm import tqdm
import polars as pl
import warnings

from subprocess import check_output

PROPERTIES_TO_EXTRACT_FROM_MOLS = ['PUBCHEM_COMPOUND_CID', 'PUBCHEM_IUPAC_NAME', 'PUBCHEM_SMILES']


def get_gzipped_sdf_filenames(directory):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('sdf.gz')]


def get_number_of_lines(filename):
    """Get the number of lines in a file. Around 2x faster than reading in all lines and counting them.

    Args:
        filename (str): The file to count the lines of

    Returns:
        int: The number of lines in the file
    """
    return int(check_output(["wc", "-l", filename]).split()[0])


def process_gzipped_sdf_file(gzipped_sdf_filename, silence_tqdm=False):
    data = {property: [] for property in PROPERTIES_TO_EXTRACT_FROM_MOLS}

    num_lines = get_number_of_lines(gzipped_sdf_filename)

    with gzip.open(gzipped_sdf_filename, 'rt') as f:
        
        next_field_to_add_to = None
        
        for line in tqdm(f, desc=f'Processing {gzipped_sdf_filename}', disable=silence_tqdm, total=num_lines):
            for prop in PROPERTIES_TO_EXTRACT_FROM_MOLS:
                if prop in line:
                    next_field_to_add_to = prop
                    break
            
            if next_field_to_add_to is not None:
                data[next_field_to_add_to].append(line.strip())
                next_field_to_add_to = None
        
    return data


def main(args):
    breakpoint()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if args.num_processes != 1:
        warnings.warn('multiprocessing is not implemented yet')

    gzipped_sdf_filenames = get_gzipped_sdf_filenames(args.input_dir)

    data_0 = process_gzipped_sdf_file(gzipped_sdf_filenames[0], silence_tqdm=False)

    data_0 = pl.from_dict(data_0)

    data_0.write_csv(os.path.join(args.output_dir, 'test.csv'))


def parse_args():
    parser = argparse.ArgumentParser(description='Converts .sdf files to .smi files')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing pubchem data as .sdf.gz')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save data')
    parser.add_argument('--num_processes', type=int, default=-1, help='Number of processes to use. -1 means use all available cores')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)