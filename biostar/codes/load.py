"""Used as an independent module to load into codes_db from command line."""

from django.core.exceptions import AppRegistryNotReady, ImproperlyConfigured

try:
    # Try import models without setting up django
    from biostar.codes.auth import load_qpcr
except (AppRegistryNotReady, ImproperlyConfigured):
    import django
    import os
    # Set up django before loading models ( only run once ).
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "biostar.settings")
    django.setup()
    from biostar.codes.auth import load_qpcr


def main():
    from argparse import ArgumentParser

    parse = ArgumentParser()

    parse.add_argument('-sample_sheet', dest='sample_sheet', required=True,
                       help="""CSV file used to load database with metadata information about samples.""",
                       type=str)
    parse.add_argument('-data_dir', dest='data_dir',
                       help="""Data directory with csv files outputted by qPCR instrument.""",
                       type=str)
    parse.add_argument('-zipfile', dest='zipfile', type=str,
                       help="""Zip filled with csv file from qPCR instrument. 
                                Unpacked to created experiments and measurement.""")

    args = parse.parse_args()
    load_qpcr(sample_sheet=args.sample_sheet, data_dir=args.data_dir, zip_file=args.zipfile)


if __name__ == "__main__":

    import sys
    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()
