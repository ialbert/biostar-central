#!/usr/bin/env python
# coding: utf-8


import requests
import os
import time
from urllib.request import urljoin
from datetime import datetime
from zipfile import ZipFile
from io import StringIO
import argparse


# Api endpoints
#ROOT_URL = "https://www.bioinformatics.recipes"

ROOT_URL = "http://localhost:8000"

PROJECT_URL = urljoin(ROOT_URL, "/api/project/")
RECIPES_URL = urljoin(ROOT_URL, "/api/recipe/")
LISTING_URL = urljoin(ROOT_URL, "/api/list/")


def show_progress():
    verbose = os.getenv('VERBOSE', False)

    return verbose


def compress(datatuple, outname=None):
    """
    Compress a list of files
    datatuple has the structure : [(fname, stream) , .... ]
    """

    stream = ZipFile(outname, 'w')
    for fname, data in datatuple:
        stream.writestr(fname, data.getvalue())
    return


def send_request(url, method='GET', params={}, data={}):
    """
    Send a request with given method and data
    """

    if method == 'GET':
        response = requests.get(url, params=params)
    else:
        response = requests.post(url, params=params, data=data)

    data = response.content.decode()
    if response.status_code == 200:
        return data
    else:
        print(f'Error:Status_Code={response.status_code},Content:{data}')

    return


def listing(endpoint, token=''):
    """
    Send GET request to listing endpoint and return a unique set of ids.
    """
    params = {'token': token}
    data = send_request(url=endpoint, params=params)
    data = data.split('\n')
    uids = set()
    for line in data:
        # Parse the first item in the line.
        uids.add(line.split('\t')[0])

    return uids


def get_data(endpoint, *uids, outdir='', token='', sleep=0):

    outputs = []
    for uid in uids:
        if show_progress():
            print(f'*** Fetching: {uid}')

        # Prepare data
        params = {'token': token, 'uid': uid}
        filename = f'{uid}.json'

        # Process request and return payload
        data = send_request(url=endpoint, params=params)
        if data:
            # Store payload
            stream = StringIO(data)
            outputs.append((filename, stream))

        # Grace period before sending the next request.
        time.sleep(sleep)

    # Construct final filename.
    now = datetime.now()
    outname = f'{now.year}-{now.month}-{now.day}.zip'
    outname = os.path.abspath(os.path.join(outdir, outname))

    # Compress .json files into one zip.
    compress(datatuple=outputs, outname=outname)

    return outname


def get_projects(token='', outdir=''):
    """
    Get all projects a user has access to.
    """
    # Get project uids
    uids = listing(endpoint=LISTING_URL, token=token)

    # Dump data
    zf = get_data(PROJECT_URL, *uids, token=token, outdir=outdir)
    return zf


# def get_recipes(token='', outdir=''):
#     """
#     Get all recipes a user has access to.
#     """
#     # Get project uids
#     uids = listing(endpoint=LISTING_URL, token=token)
#
#     # Dump data
#     zf = get_data(RECIPES_URL, *uids, token=token, outdir=outdir)
#     return zf


def main():

    parser = argparse.ArgumentParser(description='Get all your projects')
    parser.add_argument('-t', '--token', type=str,
                        help='Profile token used to access API.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress.')
    parser.add_argument('-o', '--outdir', type=str,
                        help='Output directory of the dump file', default='')

    args = parser.parse_args()

    # Get it from your profile
    token = args.token
    outdir = args.outdir

    zf = get_projects(token=token, outdir=outdir)

    print(f'*** Created: {zf}')


if __name__ == '__main__':
    main()
