"""
Filter out PubChem compounds with low Tanimoto similarity scores
"""

import os
import argparse
import numpy as np
import polars as pl
import multiprocessing as mp
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from neattime import neattime

NUM_FILES_TO_TEST = 5
NUM_MOLS_TO_TEST = 100

from extract_data_from_pubchem_sdf import PROPERTIES_TO_EXTRACT_FROM_MOLS

def filter_tanimoto_results(tanimoto_directory: str, threshold: float, test: bool = False) -> pl.DataFrame:
    """Filter PubChem compounds to only include those with at least one Tanimoto score above the specified threshold. 

    Args:
        tanimoto_directory (str): the directory containing the Tanimoto similarity results as .zst files
        threshold (float): the threshold for filtering out low scores.
        test (bool, optional): whether to run in test mode. Defaults to False.

    Returns:
        pl.DataFrame: a polars DataFrame containing the filtered Tanimoto results
    """
    tanimoto_result_filepaths = [os.path.join(tanimoto_directory, filename) for filename in os.listdir(tanimoto_directory) if filename.endswith('.zst')]

    if test:
        tanimoto_result_filepaths = tanimoto_result_filepaths[:NUM_FILES_TO_TEST]

    combine_and_filter_query = (
        pl.scan_parquet(tanimoto_result_filepaths)
        .filter(pl.any_horizontal(pl.exclude(PROPERTIES_TO_EXTRACT_FROM_MOLS) > threshold))
    )

    filtered_tanimoto_results = combine_and_filter_query.collect()

    return filtered_tanimoto_results


def main(args):
    os.makedirs(args.output_directory, exist_ok=True)

    if args.test:
        args.output_filename = f'TEST_{args.output_filename}'
        
    filtered_tanimoto_results = filter_tanimoto_results(tanimoto_directory=args.tanimoto_directory, threshold=args.threshold, test=args.test)

    filtered_tanimoto_results_save_path = os.path.join(args.output_directory, args.output_filename)

    filtered_tanimoto_results.write_csv(filtered_tanimoto_results_save_path)

    print(f'Filtered Tanimoto results saved to {filtered_tanimoto_results_save_path}')


def parse_args():
    parser = argparse.ArgumentParser(description='Filter results of tanimoto similarity search to only include PubChem compounds with at least one score above a specified threshold')
    parser.add_argument('--tanimoto_directory', type=str, required=True, help='Directory containing the Tanimoto result files')
    parser.add_argument('--output_directory', type=str, default='filtered_tanimoto_similarity_results', help='Directory to save the filtered results. Will be stored as a single csv')
    parser.add_argument('--output_filename', type=str, default=f'filtered_tanimoto_results_{neattime()}.csv', help='Filename for the filtered results')
    parser.add_argument('--threshold', type=float, default=0.8, help='Threshold for filtering out low scores. Final output will only contain scores above this threshold.')
    parser.add_argument('--test', action='store_true', help='Run in test mode (process only a subset of files and compounds)')
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    main(args)
