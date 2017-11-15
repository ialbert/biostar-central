from biostar.tools.render import render_template
import glob
import os


if __name__ == "__main__":

    iobio = "http://bam.iobio.io/?bam="

    import argparse

    parser = argparse.ArgumentParser(description='Create iobio links for bam files.')

    parser.add_argument('--base', help='The URL to file location base')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')

    args = parser.parse_args()

    # Read the arguments
    base = args.base
    patt = args.bams

    patt = os.path.join(patt, '*.bam')

    bams, links  = [], []

    for fname in glob.glob(patt):
        # print(fname)

        name = os.path.basename(fname)
        path = os.path.join(base, name)

        bams.append((path, name))

    # Create iobio link.
    for path, name in bams:
        path = path.replace(":","%3A")
        path = path.replace("/", "%2F")
        biopath = iobio + path
        links.append((biopath, name))

    # render
    data = dict(links=links)
    name = "iobio.html"

    html = render_template(data, name)
    print(html)




