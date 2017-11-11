from biostar.tools.render import render_template
import sys, os


def make_igv(genome, bam):

    bam_dir = os.path.dirname(bam)
    bam = os.path.basename(bam)

    igv_template = "igv.xml"
    data = dict(genome=genome, bam_dir=bam_dir, bam=bam)
    html = render_template(data, igv_template)
    print(html)


if __name__ == "__main__":
    genome = sys.argv[1]
    bam = sys.argv[2]
    #genome="/Users/asebastian/work/web-dev/biostar-engine/export/local/input/viral_genomes.fa"
    #bam="/Users/asebastian/work/web-dev/biostar-engine/export/local/bam/virus_sample1.bam"

    make_igv(genome, bam)


