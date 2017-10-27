from biostar.tools.plotify import ChartParams
from biostar.tools.render import render_template
import sys, csv , re


def parse_flagstats(fname):
    data





def make_plot():

    file1 = sys.argv[1]
    #file1 = "chrom_mapping.txt"
    data1, data2 = [], []

    # Plot samtools idxstats results.
    with open(file1) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            chrom = row['Chrom']

            patt = r'gi\|[0-9]{1,}\|[a-z]{1,}\|'
            chrom = re.sub(patt, '', chrom)
            chrom = chrom.replace('|', '')

            mapped = row['Mapped']
            data1.append((chrom, mapped))

    data1.pop()

    p1 = ChartParams()
    p1.type = 'BarChart'
    p1.data = data1

    p1.xlabel = "Chromosomes"
    p1.ylabel = "Mapped reads"

    p1.options = '''    
            title: 'Mapped counts per chromosome.',
            legend: {position: 'none'},
        '''




    '''
    p2 = ChartParams()
    p2.type = 'PieChart'
    p2.data = data2

    p2.xlabel = "Category"
    p2.ylabel = "Percent"


    '''
    # This is the context.
    data = dict(p1=p1)

    name = "bwa.html"

    html = render_template(data, name)
    print(html)

    # with open('index.html', 'wt') as fp:
    #    fp.write(html)


if __name__ == '__main__':
    parse_flagstats("alignment_stats.txt")
    make_plot()

