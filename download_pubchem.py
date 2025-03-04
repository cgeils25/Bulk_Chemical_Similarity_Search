import os
import requests
import time
import warnings
from neattime import neattime
import re

test = True

data_save_dir = f'pubchem_data/{'TEST_' if test else 'full_download'}{neattime()}/'

if not os.path.exists(data_save_dir):
    os.makedirs(data_save_dir)

print(f"Data will be saved to {data_save_dir}")

url_root = "https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/CURRENT-Full/SDF/"

def get_all_gzip_urls(url_root):
    response = requests.get(url_root)

    response_text = response.text

    pattern = r'<a href="(Substance_\d+_\d+\.sdf\.gz)">'
    urls = re.findall(pattern, response_text)
    full_urls = [url_root + u for u in urls]
    return full_urls

def main():
    print("Retrieving all gzip urls")
    urls = get_all_gzip_urls(url_root)

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

        if test and i > 3:
            break
    
    print("All done.")

    success_count = status_codes.count(200)
    print(f"{success_count}/{num_urls} requests were successful.")

    if success_count < num_urls:
        print('Failed requests:')
        for i, (url, status_code) in enumerate(zip(urls, status_codes)):
            if status_code != 200:
                print(f'{i}: url: {url}, status code: {status_code}')
    else:
        print('No failed requests.')

if __name__ == '__main__':
    main()
