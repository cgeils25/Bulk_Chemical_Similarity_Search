"""
--what does this file do--

To be clear, you can just download the .smi files directly from PubChem, but I wanted to have all the data in one place in case I need it for the future
"""

import os
import gzip
import argparse
import multiprocessing as mp
import polars as pl
import warnings
from tqdm import tqdm
from neattime import neattime
from rdkit import Chem

from subprocess import check_output

PROPERTIES_TO_EXTRACT_FROM_MOLS = ['PUBCHEM_COMPOUND_CID', 'PUBCHEM_IUPAC_NAME', 'PUBCHEM_SMILES']


def get_gzipped_sdf_filenames(directory):
    return [filename for filename in os.listdir(directory) if filename.endswith('sdf.gz')]


def process_gzipped_sdf_file(gzipped_sdf_filepath, test=False):

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


def save_data_to_csv(data, output_dir, filename):
    df = pl.from_dict(data)

    df.write_csv(os.path.join(output_dir, filename))

    print(f'Saved data to {os.path.join(output_dir, filename)}')


def process_gzipped_sdf_file_and_save(gzipped_sdf_filepath, output_dir, test=False):
    data = process_gzipped_sdf_file(gzipped_sdf_filepath, test=test)

    filename = os.path.basename(gzipped_sdf_filepath).replace('.sdf.gz', '.csv')

    save_data_to_csv(data, output_dir, filename)

    print(f'Processed {gzipped_sdf_filepath}')


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

            out_filename = gzipped_sdf_filename.replace('.sdf.gz', '.csv')

            save_data_to_csv(data, args.output_dir, out_filename)

            if args.test and i == 4:
                break
    
    else:
        if args.num_processes == -1:
            num_processes = mp.cpu_count()
        else:
            num_processes = args.num_processes

        print(f'Using {num_processes} processes')
        
        pool = mp.Pool(processes=num_processes)

        gzipped_sdf_filepaths = [os.path.join(args.input_dir, filename) for filename in gzipped_sdf_filenames]

        num_files = len(gzipped_sdf_filepaths)

        pool.starmap(
            process_gzipped_sdf_file_and_save,
            list(zip(gzipped_sdf_filepaths,
            [args.output_dir] * num_files,
            [args.test] * num_files))
        )


def parse_args():
    parser = argparse.ArgumentParser(description='Converts .sdf files to .smi files')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing pubchem data as .sdf.gz')
    parser.add_argument('--output_dir', type=str, required=False, default=None, help='Directory to save the extracted data')
    parser.add_argument('--num_processes', type=int, default=-1, help='Number of processes to use. -1 means use all available cores')
    parser.add_argument('--test', action='store_true', help='Test mode.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)