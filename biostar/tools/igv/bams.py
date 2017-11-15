import glob
import os
from biostar.tools import render


def create_urls(base_url, patt):
    store = []

    for fname in glob.glob(patt):

        name = os.path.basename(fname)
        path = os.path.join(base_url, name)

        store.append((path, name))

    return store


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Create IGV tracks from bam files.')

    parser.add_argument('--bamURL', help='The URL to bam file location')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')
    parser.add_argument('--bigwigURL', help='The URL to bigwig file location')
    parser.add_argument('--bigwigs', help='The path to the directory that contain bigwig files')
    parser.add_argument('--genome', help='The path to the genome')

    args = parser.parse_args()

    # Read the arguments
    bam_url = args.bamURL
    bam_dir = args.bams
    bwig_url = args.bigwigURL
    bwig_dir =args.bigwigs
    genome = args.genome

    bam_patt = os.path.join(bam_dir, '*.bam')
    bams = create_urls(bam_url, bam_patt)

    bwig_patt = os.path.join(bwig_dir, '*.bw')
    bigwigs = create_urls(bwig_url, bwig_patt)

    """
    for fname in glob.glob(patt):
        #print(fname)

        name = os.path.basename(fname)
        path = os.path.join(bam_url, name)

        bams.append((path, name))
    """

    name = "igv/igv.xml"

    data = dict(bams=bams, bigwigs=bigwigs, genome=genome)

    html = render.render_template(data, name)

    print(html)
