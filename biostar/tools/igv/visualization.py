from biostar.tools.render import render_template
import sys, os


def make_igv(genome):
    # collect all bam files.


    bam_dir = os.path.dirname(bam)
    bam = os.path.basename(bam)

    igv_template = "igv.xml"
    data = dict(genome=genome, bam_dir=bam_dir, bam=bam)
    html = render_template(data, igv_template)
    print(html)


if __name__ == "__main__":
    genome = sys.argv[1]
    make_igv(genome)