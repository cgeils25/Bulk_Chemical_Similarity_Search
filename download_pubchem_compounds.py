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

def main(args):
    data_save_dir = f'pubchem_data/{'TEST_' if args.test else 'full_download_'}{neattime()}/'

    if not os.path.exists(data_save_dir):
        os.makedirs(data_save_dir)

    print(f"Data will be saved to {data_save_dir}")

    if args.test:
        print('-'*100, f'\nRunning in test mode. Only downloading {NUM_FILES_TO_TEST} files.')
        print('-'*100)

    print("Retrieving all gzip urls")
    urls = get_all_gzip_urls(url_root)

    if args.test:
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

        file_save_path = f"{data_save_dir}{url.split('/')[-1]}"

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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_true', help=f'Run in test mode. Only download {NUM_FILES_TO_TEST} files.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(args)
