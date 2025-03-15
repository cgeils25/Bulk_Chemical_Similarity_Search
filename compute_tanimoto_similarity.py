"""
Compute Tanimoto similarity between a comparison dataset and PubChem compounds

for more info, run:
    python compute_tanimoto_similarity.py --help
"""

import os
import sys
import argparse
import numpy as np
import polars as pl
import multiprocessing as mp
import warnings
from rdkit.Chem.SaltRemover import SaltRemover 

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog('rdApp.*') # ignore rdkit's warnings

from neattime import neattime

from utils import print_args

NUM_FILES_TO_TEST = 5
NUM_MOLS_TO_TEST = 100

# will include these in the output tanimoto results
from extract_data_from_pubchem_sdf import PROPERTIES_TO_EXTRACT_FROM_MOLS

# ignores a warning that pops up when the script ends which halts everything: source_tracker: There appear to be 2 leaked semaphore objects to clean up at shutdown: {'/mp-7hozeleu', '/mp-_ni56uyz'}
#  warnings.warn(
os.environ['PYTHONWARNINGS'] = 'ignore:resource_tracker'


def smiles_list_to_fingerprint_matrix(smiles_list: list, fingerprint_size: int, radius: int, remove_salts: bool) -> np.ndarray:
    """Generate a matrix conaining morgan fingerprints from a list of SMILES strings

    Each row in the output corresponds to a fingerprint from one SMILES string

    Args:
        smiles_list (list): a list of SMILES strings
        fingerprint_size (int, optional): the length of the morgan fingerprint.
        radius (int, optional): radius of the morgan fingerprint. This determines the number of bonds away from a central atom to consider when extracting substructures to embed.
        remove_salts (bool, optional): whether or not to remove salts from the molecules.

    Returns:
        np.ndarray: a matrix of shape (len(smiles_list), fingerprint_size) containing the morgan fingerprints
    """
    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fingerprint_size)
    
    num_mols = len(smiles_list)

    fingerprint_matrix = np.zeros((num_mols, fingerprint_size), dtype=np.int8)

    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    if remove_salts:
        print('remove_salts set to True. Removing salts from molecules...')
        salt_remover = SaltRemover()
        mols = [salt_remover.StripMol(mol) for mol in mols]
        print('Salts removed from molecules.')
    
    print('Computing morgan fingerprints...')
    
    for i, mol in enumerate(mols):
        
        if mol is None:
            print(f'Warning: Invalid SMILES at index {i}: {smiles_list[i]}')
        
        else: 
            fingerprint = fingerprint_generator.GetFingerprintAsNumPy(mol)

            fingerprint_matrix[i] = np.array(fingerprint)

    return fingerprint_matrix


def compute_tanimoto_similarity_matrix(fingerprint_matrix_1: np.ndarray, fingerprint_matrix_2: np.ndarray) -> np.ndarray:
    """Compute the pairwise Tanimoto similarity between two sets of morgan fingerprints

    Index i,j of the output matrix corresponds to the Tanimoto similarity between the i-th fingerprint in fingerprint_matrix_1 and the j-th fingerprint in fingerprint_matrix_2

    Args:
        fingerprint_matrix_1 (np.ndarray): the first matrix of morgan fingerprints
        fingerprint_matrix_2 (np.ndarray): the second matrix of morgan fingerprints

    Returns:
        np.ndarray: a matrix of shape (len(fingerprint_matrix_1), len(fingerprint_matrix_2)) containing the Tanimoto similarity
    """
    # scipy's implementation of this is weird (see Jaccard and Rogers-Tanimoto). Chose to implement it myself
    C = fingerprint_matrix_1 @ fingerprint_matrix_2.T

    A = np.sum(fingerprint_matrix_1, axis=1).reshape(-1, 1)
    B = np.sum(fingerprint_matrix_2, axis=1).reshape(-1, 1)

    distance_matrix = C / (A + B.T - C)

    return distance_matrix


def get_comparison_smiles_list(filepath: str) -> list:
    """Get the list of SMILES strings from a comparison dataset

    Args:
        filepath (str): the path to the comparison dataset. Must be a CSV file with a 'smiles' column

    Raises:
        ValueError: if the file is not a CSV file
        ValueError: if the column 'smiles' is not found in the dataset

    Returns:
        list: a list of SMILES strings
    """
    if not filepath.endswith('.csv'):
        raise ValueError(f'File {filepath} is not a CSV file')
    
    df = pl.read_csv(filepath)

    if 'smiles' not in df.columns:
        raise ValueError(f'Column "smiles" not found in {filepath}')

    smiles_list = df['smiles'].to_list()

    return smiles_list


def extracted_pubchem_data_filename_to_tanimoto_similarity_filename(extracted_pubchem_data_filename) -> str:
    extracted_pubchem_data_filename = extracted_pubchem_data_filename.replace('.zst', '_tanimoto_similarity.zst')

    return extracted_pubchem_data_filename

def get_incomplete_tanimoto_similarity_files(extracted_pubchem_data_dir: str, output_dir: str) -> list:
    print(f'Checking for incomplete tanimoto similarity files in {output_dir} based on files found in {extracted_pubchem_data_dir}...')

    extracted_pubchem_data_filenames = os.listdir(extracted_pubchem_data_dir)

    tanimtoto_similarity_filenames = os.listdir(output_dir)

    incomplete_tanimoto_similarity_files = []

    for extracted_pubchem_data_filename in extracted_pubchem_data_filenames:
        expected_tanimoto_similarity_filename = extracted_pubchem_data_filename_to_tanimoto_similarity_filename(extracted_pubchem_data_filename=extracted_pubchem_data_filename)

        if expected_tanimoto_similarity_filename not in tanimtoto_similarity_filenames:
            print(f'Found incomplete tanimoto similarity file for {extracted_pubchem_data_filename}. Expected {expected_tanimoto_similarity_filename} but not found in {output_dir}. Adding to list of files to process.')
            extracted_pubchem_data_filepath = os.path.join(extracted_pubchem_data_dir, extracted_pubchem_data_filename)

            incomplete_tanimoto_similarity_files.append(extracted_pubchem_data_filepath)
        
        else:
            print(f'Found completed tanimoto similarity file for {extracted_pubchem_data_filename} at {expected_tanimoto_similarity_filename}. Skipping.')
    
    return incomplete_tanimoto_similarity_files


def run_comparison(comparison_smiles_list: list, comparison_dataset: str, extracted_pubchem_data_filepath: str, 
                   output_dir: str, fingerprint_size: int, radius: int, remove_salts: bool, test: bool = False) -> None:
        """Run a Tanimoto similarity search between a comparison dataset and PubChem compounds and save the results

        Args:
            comparison_smiles_list (list): a list of SMILES strings from the dataset to compare PubChem against
            comparison_dataset (str): path to the comparison dataset
            extracted_pubchem_data_filepath (str): path to the extracted PubChem data
            output_dir (str): path to the directory where the results will be saved
            fingerprint_size (int): size of the morgan fingerprint
            radius (int): radius of the morgan fingerprint
            remove_salts (bool): whether or not to remove salts from the molecules
            test (bool, optional): whether to run in test mode (will only process a small sample of the data). Defaults to False.
        """
        print(f'Running tanimoto similarity computation between {extracted_pubchem_data_filepath} and {comparison_dataset}')
        extracted_pubchem_data_df = pl.read_parquet(source=extracted_pubchem_data_filepath)
        
        pubchem_smiles_list = extracted_pubchem_data_df['PUBCHEM_SMILES'].to_list()
        
        if test:
            pubchem_smiles_list = pubchem_smiles_list[:NUM_MOLS_TO_TEST]
        
        pubchem_fingerprint_matrix = smiles_list_to_fingerprint_matrix(smiles_list=pubchem_smiles_list, fingerprint_size=fingerprint_size, 
                                                                       radius=radius, remove_salts=remove_salts)
        
        comparison_fingerprint_matrix = smiles_list_to_fingerprint_matrix(smiles_list=comparison_smiles_list, fingerprint_size=fingerprint_size,
                                                                          radius=radius, remove_salts=remove_salts) # This gets computed multiple times because I don't think you can pickle it for multiprocessing. Idk ig I should test at some point
        
        tanimoto_similarity_matrix = compute_tanimoto_similarity_matrix(pubchem_fingerprint_matrix, comparison_fingerprint_matrix)
        
        tanimoto_similarity_df = pl.DataFrame(tanimoto_similarity_matrix, schema = comparison_smiles_list)
        # add the properties that were extracted from the original pubchem sdf files
        
        for i, property in enumerate(PROPERTIES_TO_EXTRACT_FROM_MOLS):
            if property in extracted_pubchem_data_df.columns:
                property_values = extracted_pubchem_data_df[property].to_list()

                if test:
                    property_values = property_values[:NUM_MOLS_TO_TEST]

                tanimoto_similarity_df.insert_column(index=i, column=pl.Series(name=property, values=property_values))
            else:
                print(f'Warning: Property {property} not found in {extracted_pubchem_data_filepath}. Skipping.')

        extracted_pubchem_data_filename = os.path.basename(extracted_pubchem_data_filepath)

        tanimoto_similarity_filename = extracted_pubchem_data_filename_to_tanimoto_similarity_filename(extracted_pubchem_data_filename=extracted_pubchem_data_filename)

        tanimoto_similarity_save_path = os.path.join(output_dir, tanimoto_similarity_filename)

        # save a compressed parquet file. My experiments show about a 2x compression rate for the test run with zstd and compression level 8 
        tanimoto_similarity_df.write_parquet(file=tanimoto_similarity_save_path, compression = 'zstd', compression_level=8) # zstd compression recommended by polars for good compression performance: https://docs.pola.rs/api/python/dev/reference/api/polars.DataFrame.write_parquet.html

        print(f'Saved tanimoto similarity between {extracted_pubchem_data_filepath} and {comparison_dataset} to {tanimoto_similarity_save_path}')
        

def main(args):
    if args.test:
        print('-'*100)
        print(f'Running in test mode. Only processing {NUM_FILES_TO_TEST} files and {NUM_MOLS_TO_TEST} molecules')
    
    print('-'*100)
    print('Running Tanimoto similarity search')
    
    print_args(args)

    comparison_smiles_list = get_comparison_smiles_list(args.comparison_dataset)

    if args.test:
        comparison_smiles_list = comparison_smiles_list[:NUM_MOLS_TO_TEST]

    if args.output_dir is None:
        args.output_dir = f'tanimoto_similarity_results/{'TEST_' if args.test else 'full_tanimoto_'}{neattime()}'

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # get the ones that aren't already done
    if args.cleanup:
        extracted_pubchem_data_filepaths = get_incomplete_tanimoto_similarity_files(extracted_pubchem_data_dir=args.extracted_pubchem_data_dir, output_dir=args.output_dir)

        if len(extracted_pubchem_data_filepaths) == 0:

            print(f'No incomplete tanimoto similarity files found in {args.output_dir}. Exiting.')

            sys.exit(0)
    
    # if not cleanup, just prepare to compute similarity for all of them
    else:
        extracted_pubchem_data_filenames = os.listdir(args.extracted_pubchem_data_dir)

        extracted_pubchem_data_filepaths = [os.path.join(args.extracted_pubchem_data_dir, filename) for filename in extracted_pubchem_data_filenames]

    if args.num_processes == 1:
        print('Using a single process')

        for i, extracted_pubchem_data_filepath in enumerate(extracted_pubchem_data_filepaths):

            run_comparison(comparison_smiles_list=comparison_smiles_list, comparison_dataset=args.comparison_dataset, 
                           extracted_pubchem_data_filepath=extracted_pubchem_data_filepath, output_dir=args.output_dir, fingerprint_size=args.fingerprint_size, 
                           radius=args.radius, remove_salts=args.remove_salts, test=args.test)
            
            if args.test and i == NUM_FILES_TO_TEST - 1:
                break

    else:
        processes_available = mp.cpu_count()
        if args.num_processes == -1:
            num_processes = processes_available
        else:
            if args.num_processes > processes_available:
                warnings.warn(f'Number of processes {args.num_processes} is greater than the number of available processes {processes_available}. Using {processes_available} instead.')
            
            num_processes = args.num_processes

        print(f'Using {num_processes} processes')

        # is there a less dumb way to do this?
        num_pubchem_files = len(extracted_pubchem_data_filepaths)

        inputs_zipped = list(zip([comparison_smiles_list] * num_pubchem_files,
                                [args.comparison_dataset] * num_pubchem_files,
                                 extracted_pubchem_data_filepaths, 
                                 [args.output_dir] * num_pubchem_files, 
                                 [args.fingerprint_size] * num_pubchem_files, 
                                 [args.radius] * num_pubchem_files, 
                                [args.remove_salts] * num_pubchem_files,
                                 [args.test] * num_pubchem_files))
        
        if args.test:
            inputs_zipped = inputs_zipped[:NUM_FILES_TO_TEST]
        
        # have to use spawn. Polars breaks if fork is used: https://docs.pola.rs/user-guide/misc/multiprocessing/#summary
        with mp.get_context("spawn").Pool(processes=num_processes) as pool:

            pool.starmap(run_comparison, inputs_zipped)
        

def parse_args():
    parser = argparse.ArgumentParser(description='''Compute pairwise tanitotomo similarity between a comparison dataset and PubChem compounds, then save results as multiple .zst files
                                     
                                     NOTE: be careful with the number of processes for `compute_tanimoto_similarity.py`. When running this on an HPC node with a 48-core intel xeon and 196 GB of RAM, the highest I could do while remaining stable was 16. 
                                     This is the most memory-intense part of the pipeline so adjust to your needs. It seems that polars becomes unstable when combined with python multiprocessing.
                                     ''')
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
        help="Path to the directory containing PubChem data extracted from sdf files. Expected to contain multiple zst files with a PUBCHEM_SMILES column."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default=None,
        help="Path to the output directory where the similarity data will be saved. If not specified, will be saved to tanimoto_similarity_results/ with a timestamp.",
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
        "--remove_salts",
        action="store_true",
        help="Whether to remove salts from BOTH the PubChem molecules and the comparison dataset molecules prior to generating morgan fingerprints and computing tanimoto similarity. This is recommended if you are searching for similar organics.",
    )
    parser.add_argument(
        '--num_processes', 
        type=int, 
        default=-1, 
        help='Number of processes to use. -1 means use all available cores'
    )
    parser.add_argument(
        '--cleanup',
        action='store_true',
        help='Whether to search the output directory for completed tanimoto similarity results and only compute the ones that are missing. Useful for resuming a job that was interrupted or froze.'
    )
    parser.add_argument(
        '--test', 
        action='store_true', 
        help=f'Test mode. Only process {NUM_FILES_TO_TEST} files and {NUM_MOLS_TO_TEST} molecules'
    )
    
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

    