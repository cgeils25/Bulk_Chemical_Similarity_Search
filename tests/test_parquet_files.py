"""
Tests the integrity of parquet files in a specified directory by attempting to read each file.
"""

import os
import argparse
import polars as pl


def main(dirname: str):
    filenames = os.listdir(dirname)

    filepaths = [os.path.join(dirname, filename) for filename in filenames]

    successes = 0
    failures = 0

    for filepath in filepaths:
        try:
            # Read the parquet file
            pl.read_parquet(filepath)
            print(f"Successfully read {filepath}")
            successes += 1
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            failures += 1

    print(f"Successfully read {successes} files.")
    print(f"Failed to read {failures} files.")

def parse_args():
    parser = argparse.ArgumentParser(description='Test the integrity of parquet files in a specified directory by attempting to read each file.')
    parser.add_argument('dirname', type=str, help='The directory containing the parquet files to test.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(dirname = args.dirname)
    
