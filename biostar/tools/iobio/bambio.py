from biostar.tools.render import render_template
import os, glob


def bam_links(base_url, bam_dir):

    patt = os.path.join(bam_dir, '*.bam')

    bams = []

    for fname in glob.glob(patt):

        name = os.path.basename(fname)
        path = os.path.join(base_url, name)
        path = path.replace(":", "%3A")
        path = path.replace("/", "%2F")

        bams.append((path, name))

    return bams


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Create iobio links for bam files.')

    parser.add_argument('--base', help='The URL to file location base')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')

    args = parser.parse_args()

    # Read the arguments.
    base = args.base
    bamdir = args.bams

    # Create iobio bam links.
    bams = bam_links(base, bamdir)

    # Render template.
    data = dict(bams=bams)
    name = "align/scaffold_index.html"

    html = render_template(data, name)
    print(html)




