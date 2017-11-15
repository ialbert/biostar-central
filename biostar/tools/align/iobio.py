from biostar.tools.render import render_template
from biostar.tools.igv.paths import collect_paths


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Create iobio links for bam files.')

    parser.add_argument('--base', help='The URL to file location base')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')

    args = parser.parse_args()

    # Read the arguments
    base = args.base
    bamdir = args.bams
    patt = "*.bam"

    bams = collect_paths(base, bamdir, patt)

    links = []
    iobio = "http://bam.iobio.io/?bam="

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




