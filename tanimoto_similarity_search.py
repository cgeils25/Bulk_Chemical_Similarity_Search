"""
Run a Tanimoto similarity search
"""

import os
import argparse
import numpy as np
import polars as pl
import multiprocessing as mp
from tqdm import tqdm
import gzip

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from neattime import neattime

NUM_FILES_TO_TEST = 5
NUM_MOLS_TO_TEST = 100


def smiles_list_to_fingerprint_matrix(smiles_list, fingerprint_size=2048, radius=2, disable_tqdm=False):
    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fingerprint_size)
    
    num_mols = len(smiles_list)

    fingerprint_matrix = np.zeros((num_mols, fingerprint_size), dtype=np.int8)

    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    
    for i, mol in enumerate(tqdm(mols, desc='Computing fingerprints', disable=disable_tqdm)):
        
        if mol is None:
            print(f'Warning: Invalid SMILES at index {i}: {smiles_list[i]}')
        
        else: 
            fingerprint = fingerprint_generator.GetFingerprintAsNumPy(mol)

            fingerprint_matrix[i] = np.array(fingerprint)

    return fingerprint_matrix


def compute_tanimoto_similarity_matrix(fingerprint_matrix_1, fingerprint_matrix_2):
    # scipy's implementation of this is weird (see Jaccard and Rogers-Tanimoto). Chose to implement it myself
    C = fingerprint_matrix_1 @ fingerprint_matrix_2.T

    A = np.sum(fingerprint_matrix_1, axis=1).reshape(-1, 1)
    B = np.sum(fingerprint_matrix_2, axis=1).reshape(-1, 1)

    distance_matrix = C / (A + B.T - C)

    return distance_matrix


def get_comparison_smiles_list(filepath):
    df = pl.read_csv(filepath)

    if 'smiles' not in df.columns:
        raise ValueError(f'Column "smiles" not found in {filepath}')

    smiles_list = df['smiles'].to_list()

    return smiles_list


def run_comparison(comparison_smiles_list, comparison_dataset, extracted_pubchem_data_filepath, output_dir, fingerprint_size, radius, test=False, disable_tqdm=False):
        extracted_pubchem_data_df = pl.read_csv(source=extracted_pubchem_data_filepath)

        pubchem_smiles_list = extracted_pubchem_data_df['PUBCHEM_SMILES'].to_list()

        if test:
            pubchem_smiles_list = pubchem_smiles_list[:NUM_MOLS_TO_TEST]

        pubchem_fingerprint_matrix = smiles_list_to_fingerprint_matrix(pubchem_smiles_list, fingerprint_size, radius, disable_tqdm=disable_tqdm)
        comparison_fingerprint_matrix = smiles_list_to_fingerprint_matrix(comparison_smiles_list, fingerprint_size, radius, disable_tqdm=disable_tqdm) # This gets computed multiple times because I don't think you can pickle it for multiprocessing. Idk ig I should test at some point

        tanimoto_similarity_matrix = compute_tanimoto_similarity_matrix(pubchem_fingerprint_matrix, comparison_fingerprint_matrix)

        tanimoto_similarity_df = pl.DataFrame(tanimoto_similarity_matrix, schema = comparison_smiles_list)

        tanimoto_similarity_df.insert_column(index=0, column=pl.Series(name='PUBCHEM_SMILES', values=pubchem_smiles_list))

        extracted_pubchem_data_filename = os.path.basename(extracted_pubchem_data_filepath).replace('.csv', '_tanimoto_similarity.zst')

        save_path = os.path.join(output_dir, extracted_pubchem_data_filename) 

        # save a compressed parquet file. My experiments show about a 2x compression rate for the test run with zstd and compression level 8 
        tanimoto_similarity_df.write_parquet(file=save_path, compression = 'zstd', compression_level=8) # zstd compression recommended by polars for good compression performance: https://docs.pola.rs/api/python/dev/reference/api/polars.DataFrame.write_parquet.html

        print(f'Saved tanimoto similarity between {extracted_pubchem_data_filepath} and {comparison_dataset} to {save_path}')


def run_comparison_parallel_support(comparison_dataset, extracted_pubchem_data_filepath, output_dir, fingerprint_size, radius, test=False):
    comparsion_smiles_list = get_comparison_smiles_list(comparison_dataset)

    run_comparison(comparsion_smiles_list, comparison_dataset, extracted_pubchem_data_filepath, output_dir, fingerprint_size, radius, test, disable_tqdm=True)


def main(args):
    comparison_smiles_list = get_comparison_smiles_list(args.comparison_dataset)

    if args.test:
        comparison_smiles_list = comparison_smiles_list[:NUM_MOLS_TO_TEST]

    if args.output_dir is None:
        args.output_dir = f'tanimoto_similarity_results/{'TEST_' if args.test else 'full_tanimoto_'}{neattime()}'
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    extracted_pubchem_data_filenames = os.listdir(args.extracted_pubchem_data_dir)

    extracted_pubchem_data_filepaths = [os.path.join(args.extracted_pubchem_data_dir, filename) for filename in extracted_pubchem_data_filenames]

    if args.num_processes == 1:
        print('Using a single process')

        for i, extracted_pubchem_data_filepath in enumerate(extracted_pubchem_data_filepaths):

            run_comparison(comparison_smiles_list, args.comparison_dataset, extracted_pubchem_data_filepath, args.output_dir, args.fingerprint_size, args.radius, args.test)
            
            if args.test and i == NUM_FILES_TO_TEST - 1:
                break

    else:
        raise NotImplementedError('Multiprocessing does not work yet. Just freezes on me. Confirmed tqdm was not the problem')
        if args.num_processes == -1:
            num_processes = mp.cpu_count()
        else:
            num_processes = args.num_processes

        print(f'Using {num_processes} processes')

        # is there a less dumb way to do this?
        num_pubchem_files = len(extracted_pubchem_data_filenames)

        inputs_zipped = list(zip([args.comparison_dataset] * num_pubchem_files,
                                 extracted_pubchem_data_filepaths, 
                                 [args.output_dir] * num_pubchem_files, 
                                 [args.fingerprint_size] * num_pubchem_files, 
                                 [args.radius] * num_pubchem_files, 
                                 [args.test] * num_pubchem_files))
        
        with mp.Pool(processes=num_processes) as pool:

            pool.starmap(run_comparison_parallel_support, inputs_zipped)


def parse_args():
    parser = argparse.ArgumentParser(description="Run a Tanimoto similarity search")
    parser.add_argument(
        "--comparison_dataset",
        type=str,
        required=True,
        help="Path to the comparison dataset. Must be a CSV file with a 'smiles' column.",
    )
    parser.add_argument(
        "--extracted_pubchem_data_dir",
        type=str,
        required=True,
        help="Path to the directory containing PubChem data extracted from sdf files. Expected to contain multiple CSV files with a PUBCHEM_SMILES column."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default=None,
        help="Path to the output directory where the similarity data will be saved.",
    )
    parser.add_argument(
        "--fingerprint_size",
        type=int,
        default=2048,
        help="Size of the morgan fingerprint",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=3,
        help="Radius of the morgan fingerprint aka the number of bonds away from a central atom to consider when extracting substructures to embed",
    )
    parser.add_argument(
        '--num_processes', 
        type=int, 
        default=-1, 
        help='Number of processes to use. -1 means use all available cores'
    )
    parser.add_argument(
        '--test', 
        action='store_true', 
        help='Test mode.'
    )
    
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

    