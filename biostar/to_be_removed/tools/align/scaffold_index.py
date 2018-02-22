from biostar.tools.align.scripts import scaffold_plotter
from biostar.tools.iobio import bambio
from biostar.tools.render import render_template

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create barcharts from samtools idxstats output.')

    parser.add_argument('--base', help='The URL to bam file location base')
    parser.add_argument('--bams', help='The path to the directory that contain bam files')
    parser.add_argument('--mapped', help='File with idxstats output ')
    parser.add_argument('--total', help='File with total read counts ')
    parser.add_argument('--selected', help='No. of reads subselected for mapping')

    args = parser.parse_args()

    # Read the arguments
    bambase = args.base
    bamdir = args.bams
    mapped = args.mapped
    total = args.total
    selected = args.selected

    # Create plot.
    plot1 = scaffold_plotter.create_plot(mapped, total, selected)

    # Create iobio links.
    bams = bambio.bam_links(bambase, bamdir)

    # This is the context
    data = dict(p1=plot1, bams=bams)

    # Render template.
    name = "align/scaffold_index.html"

    html = render_template(data, name)
    print(html)



