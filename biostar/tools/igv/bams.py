import glob
import os

from biostar.tools import render

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Create IGV tracks from bam files.')

    parser.add_argument('--base', help='The URL to file location base')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')
    parser.add_argument('--genome', help='The path to the genome')

    args = parser.parse_args()

    # Read the arguments
    base = args.base
    patt = args.bams
    genome = args.genome

    patt = os.path.join(patt, '*.bam')

    bams = []

    for fname in glob.glob(patt):
        #print(fname)

        name = os.path.basename(fname)
        path = os.path.join(base, name)

        bams.append((path, name))

    name = "igv/igv.xml"

    data = dict(bams=bams, genome=genome)

    html = render.render_template(data, name)

    print(html)
