from biostar.tools.plotify import ChartParams
from biostar.tools.render import render_template
import sys, csv , re


def make_plot():

    fname = sys.argv[1]
    #fname = "chrom_mapping.txt"
    items = []

    # Read alignment report file.
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            chrom = row['Chrom']

            patt = r'gi\|[0-9]{1,}\|[a-z]{1,}\|'
            chrom = re.sub(patt, '', chrom)
            chrom = chrom.replace('|', '')

            mapped = row['Mapped']
            items.append((chrom, mapped))

    items.pop()

    p1 = ChartParams()
    p1.type = 'BarChart'
    p1.data = items

    p1.xlabel = "Chromosomes"
    p1.ylabel = "Mapped reads"

    p1.options = '''    
            title: 'Mapped counts per chromosome.',
            legend: {position: 'none'},
        '''

    # This is the context.
    data = dict(p1=p1)

    name = "alignment.html"

    html = render_template(data, name)
    print(html)

    #with open('index.html', 'wt') as fp:
    #    fp.write(html)


if __name__ == '__main__':
    make_plot()

