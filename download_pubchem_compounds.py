"""
Download all molecules from the PubChem Compound database as SDF files.

Use flag --test to only download a few files for testing purposes.
"""

import os
import requests
import time
from neattime import neattime
import re
import argparse

from utils import print_args

NUM_FILES_TO_TEST = 3

url_root = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"

def get_all_gzip_urls(url_root: str) -> list:
    """Get all urls ending with .sdf.gz from the given url_root.

    Args:
        url_root (str): url to the directory containing the .sdf.gz files

    Returns:
        list: list of urls ending with .sdf.gz
    """

    response = requests.get(url_root)

    response_text = response.text

    pattern = r'<a href="(Compound_\d+_\d+\.sdf\.gz)">'
    urls = re.findall(pattern, response_text)
    full_urls = [url_root + u for u in urls]
    return full_urls

def download_pubchem_compounds(urls: list, output_dir: str, test: bool) -> None:
    if test:
        print(f'Running in test mode. Only downloading {NUM_FILES_TO_TEST} files.')
        urls = urls[:NUM_FILES_TO_TEST]

    num_urls = len(urls)

    print(f"Extracted {num_urls} urls.")

    status_codes = []

    for i, url in enumerate(urls):
        print('-'*100, f'\nDownloading url #{i}: {url}')

        print('Making request...')
        response = requests.get(url)

        status_code = response.status_code

        print(f'Request status code: {status_code}')

        status_codes.append(status_code)

        print('Saving data...')

        file_save_path = f"{output_dir}{url.split('/')[-1]}"

        with open(file_save_path, 'wb') as f:
            f.write(response.content)

        print(f'Data saved to {file_save_path}')

        # avoid making too many requests to quickly
        time.sleep(1)
    
    print("All done.")

    success_count = status_codes.count(200)
    print(f"{success_count}/{num_urls} requests were successful.")

    if success_count < num_urls:
        print('Failed requests:')
        for i, (url, status_code) in enumerate(zip(urls, status_codes)):
            if status_code != 200:
                print(f'#{i}: url: {url}, status code: {status_code}')
    else:
        print('No failed requests.')


def main(args):
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = f'pubchem_data/{'TEST_' if args.test else 'full_download_'}{neattime()}/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Data will be saved to {output_dir}")

    if args.test:
        print('-'*100, f'\nRunning in test mode. Only downloading {NUM_FILES_TO_TEST} files.')
        print('-'*100)
    
    print_args(args)

    print("Retrieving all gzip urls")
    urls = get_all_gzip_urls(url_root)

    download_pubchem_compounds(urls, output_dir, args.test)


def parse_args():
    parser = argparse.ArgumentParser(description=f'''Download all PubChem compounds as .sdf.gz files from the pubchem ftp site ({url_root})
                                     
                                     FYI this ended up being 111 GB of compressed (.gz) files for me. Uncompressed, I calculated it would be 909.8 GB, although my implementation leaves the files themselves compressed while processing them.\n\n

                                     Add a -u flag (python -u ...) if you're writing to a log file to ensure that the output is not buffered. Otherwise, python will just print may just print everything at the script's exit
                                     ''')
    parser.add_argument('--test', action='store_true', help=f'Run in test mode. Only download {NUM_FILES_TO_TEST} files.')
    parser.add_argument('--output_dir', type=str, default=None, help='Directory to save the downloaded files. If not specified, will be saved to pubchem_data/ with a timestamp.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(args)
