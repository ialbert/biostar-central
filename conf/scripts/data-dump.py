#!/usr/bin/env python
# coding: utf-8

"""

user cases:
1. Migrating a whole project based on project ID, info and all the recipes.

input: project id
output: a file

"""
import requests
import os
import time
import sys
import argparse
from urllib.parse import urlsplit, urljoin


def resolve_url(base, view):

    mapper = dict(projects='api/project/',
                  recipes='api/recipe/',
                  listing='api/list/')

    endpoint = mapper.get(view, 'api/project/')
    api_url = urljoin(base, endpoint)

    return api_url


def send_request(url, method='GET', params={}, data={}, files={}):
    """
    Send a request with given method and data
    """

    if method == 'GET':
        response = requests.get(url, params=params)
    else:
        response = requests.post(url, params=params, data=data, files=files)

    data = response.content.decode()
    if response.status_code == 200:
        return data
    else:
        print(f'*** {response.status_code} Error: {data}')

    return


def get_data(url, uid, outdir='', token='', verbose=False):
    """
    Send a GET request to 'url' with 'uid' and 'token'.
    Saves returned data into 'outdir'
    """
    if verbose:
        print(f'*** Fetching: {uid} from {url}')

    # Prepare data
    params = {'token': token, 'uid': uid}
    filename = os.path.abspath(os.path.join(outdir, f'{uid}.json'))

    # Process request and return payload
    data = send_request(url=url, params=params)
    if data:
        # TODO: Store into gzip
        open(filename, 'w').write(data)
        if verbose:
            print(f'*** Created: {filename}')

    return filename


def listing(url, token=''):
    """
    Send GET request to listing endpoint and return a unique set of ids.
    """
    params = {'token': token}
    data = send_request(url=url, params=params)
    data = data.split('\n')
    uids = set()
    for line in data:
        # Parse the first item in the line.
        uids.add(line.split('\t')[0])

    return uids


def fetch_all(url, token, outdir='', verbose=False):
    """
    Get all projects found in remote host
    """

    # Get the base url

    # Construct the listing url
    list_api = resolve_url(base=url, view='listing')
    prj_api = resolve_url(base=url, view='projects')

    uids = listing(url=list_api, token=token)

    # Iterate over list of project ids and save to file
    for uid in uids:
        get_data(url=prj_api, uid=uid, token=token, outdir=outdir, verbose=verbose)

        # Grace period between requests
        time.sleep(3)

    return


def push_data(url, fname, uid, token=''):
    """
    Push data in --fname to object with --uid on remote host --url.
    If the --uid does not exist in remote host, a new object is created.
    """

    stream = open(fname, "r")
    files = {"data": stream}
    params = {'token': token, 'uid': uid}

    response = send_request(url=url, params=params, files=files)

    if response:
        print(f"*** Pushed {fname} to {url}")

    stream.close()

    return


def main():

    parser = argparse.ArgumentParser(description='Get all your projects')
    parser.add_argument('-t', '--token', type=str,
                        help='Profile token used to access API.')
    parser.add_argument('-u', '--url', type=str, required=True,
                        help='API endpoint.')
    parser.add_argument('--pid', type=str,
                        help='Project uid to dump.')
    parser.add_argument('--all', action='store_true', default=False,
                        help='Dump all projects in remote host.')
    parser.add_argument('--rid', type=str,
                        help='Recipe uid to dump.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress.')
    parser.add_argument('-p', '--push', action='store_true',
                        help='Push --file to the given uid.')
    parser.add_argument('-d', '--data', type=str,
                        help='File to push to remote host.')
    parser.add_argument('-o', '--outdir', type=str,
                        help='Output directory of the dump file', default='')

    args = parser.parse_args()

    # Get it from your profile
    token = args.token
    outdir = args.outdir
    url = args.url
    get_all = args.all
    verbose = args.verbose
    push = args.push
    data = args.data

    # Resolve the uid, projects always wins.
    uid = args.pid or args.rid

    # Fetch all projects from a remote host.
    if get_all:
        fetch_all(url=url, token=token, outdir=outdir, verbose=verbose)
        return

    # Resolve the endpoint using uid.
    rec_url = resolve_url(url, 'recipes')
    prj_url = resolve_url(url, 'projects')
    url = rec_url if args.rid else prj_url

    # Push argument requires --fname
    if push and not data:
        print("--data required with --push flag.")
        sys.exit()

    # Push given data to specific project/recipe  --url
    if push:
        push_data(url=url, fname=data, uid=uid, token=token)
        return

    # Get a specific project or recipe.
    get_data(url=url, uid=uid, outdir=outdir, token=token, verbose=verbose)


if __name__ == '__main__':
    main()
