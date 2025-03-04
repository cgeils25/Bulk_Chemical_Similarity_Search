import os
import requests
import time
from neattime import neattime
import re
import argparse

num_files_for_test = 3

url_root = "https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/CURRENT-Full/SDF/"

def get_all_gzip_urls(url_root):
    response = requests.get(url_root)

    response_text = response.text

    pattern = r'<a href="(Substance_\d+_\d+\.sdf\.gz)">'
    urls = re.findall(pattern, response_text)
    full_urls = [url_root + u for u in urls]
    return full_urls

def main(args):
    data_save_dir = f'pubchem_data/{'TEST_' if args.test else 'full_download_'}{neattime()}/'

    if not os.path.exists(data_save_dir):
        os.makedirs(data_save_dir)

    print(f"Data will be saved to {data_save_dir}")

    if args.test:
        print('-'*100, f'\nRunning in test mode. Only downloading {num_files_for_test} files.')
        print('-'*100)

    print("Retrieving all gzip urls")
    urls = get_all_gzip_urls(url_root)

    if args.test:
        urls = urls[:num_files_for_test]

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
    parser.add_argument('--test', action='store_true', help=f'Run in test mode. Only download {num_files_for_test} files.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(args)
