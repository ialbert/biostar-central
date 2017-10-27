from biostar.tools.plotify import ChartParams
from biostar.tools.render import render_template
import sys, csv , re


def parse_flagstats(fname):
    '''
    parse samtools flagstat results  and returns a list of (tag,value) pair.
    '''
    to_remove =["in","with"]
    patt = r'\(.*\)'
    stats = []
    stream = open(fname).readlines()
    for data in stream:
        arr = data.split('+')
        value = arr[0]
        tag = arr[1].split(" 0 ")[1]
        tag = ' '.join(filter(lambda x: x not in to_remove,  tag.split()))
        tag = re.sub(patt, "", tag)
        tag = re.sub(r'\(QC-passed reads', "", tag)
        stats.append((tag.strip(),value.strip()))
    return stats


def parse_idxstats(fname):
    '''
    parse samtools idxstats results  and returns a list of (tag,value) pair.
    '''
    stats = []
    stream = open(fname).readlines()
    for data in stream:
        arr = data.split("\t")
        chrom = arr[0]
        mapped = arr[2]

        # This is for just for fish sequences.
        patt = r'gi\|[0-9]{1,}\|[a-z]{1,}\|'
        chrom = re.sub(patt, '', chrom)
        chrom = chrom.replace('|', '')
        stats.append((chrom, mapped))

    stats.pop()
    return stats


def plot(data1, data2):

    # Sort by counts.
    data1.sort(key=lambda x: int(x[1]), reverse=True)

    # Get  only single-end specific details.
    details = ['total', 'mapped', 'supplementary', 'secondary', 'duplicates','singletons']
    data1 = list(filter(lambda x:x[0].lower() in details, data1))

    p1 = ChartParams()
    p1.type = 'BarChart'

    p1.data = data1

    p1.xlabel = "Flag category"
    p1.ylabel = "Read counts"

    p1.options = '''    
            title: 'Alignment details.',
            legend: {position: 'none'},
        '''

    p2 = ChartParams()
    p2.type = 'BarChart'

    p2.data = data2

    p2.xlabel = "Chromosomes"
    p2.ylabel = "Read counts"

    p2.options = '''    
            title: 'Mapped counts per chromosome.',
            legend: {position: 'none'},
        '''

    # This is the context.
    data = dict(p1=p1, p2=p2)

    name = "bwa.html"

    html = render_template(data, name)
    print(html)

    #with open('index.html', 'wt') as fp:
    #    fp.write(html)


if __name__ == '__main__':

    fname1 = sys.argv[1]
    fname2 = sys.argv[2]

    flag_stats = parse_flagstats(fname1)
    idx_stats = parse_idxstats(fname2)
    plot(flag_stats, idx_stats)

    #flag_stats = parse_flagstats("alignment_stats.txt")
    #idx_stats = parse_idxstats("chrom_mapping.txt")
    #plot(flag_stats, idx_stats)
