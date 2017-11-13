from biostar.tools.plotify import ChartParams,BarchartParams
from biostar.tools.render import render_template
import sys, csv , re
import csv


def get_factor():
    return 1


def parse_idxstats(fname):
    '''
    parse samtools idxstats results  and returns a list of (tag,value) pair.
    '''
    sample = None
    store = dict()
    stream = open(fname).readlines()

    for row in csv.reader(stream, delimiter='\t'):

        if len(row) == 1:
            sample = row[0]
            store[sample] = {}
            continue

        elem = store[sample]
        name, size, mapped, unmapped = row
        mapped_cov = mapped * get_factor()
        elem.setdefault(name, []).append((int(size), int(mapped_cov)))

    return store


def format_idxstats(store):
    header= ["'{0}'".format("chrom")]
    chroms = dict()

    for sample, values in store.items():
        header.append("'{0}'".format(sample) )

        for chr, stats in values.items():
            chroms.setdefault(chr, []).append(stats[0][1])

    rows = list(map(list, chroms.items()))

    # flatten list
    flattened =[]
    for elem in rows:
        arr =[]
        arr.append("'{0}'".format(elem[0]))
        for e in elem[1]:
            arr.append(e)
        flattened.append(arr)

    rows = flattened
    rows.insert(0, header)

    # remove unmapped.
    rows.pop()
    return rows


def plot(data1):

    # Plot 2
    p1 = BarchartParams()
    p1.type = 'BarChart'
    p1.header = data1[0]
    p1.data = data1[1:]
    p1.options = '''    
            title: 'Mapped counts per chromosome.',
        '''

    # This is the context.
    data = dict(p1=p1)

    name = "scaffold.html"

    html = render_template(data, name)
    print(html)


if __name__ == '__main__':

    fname1 = "mapping_stats.txt"
    #fname1 = sys.argv[1]
    #fname2 = sys.argv[2]
    store = parse_idxstats(fname1)
    idx_stats = format_idxstats(store)

    #print(store)

    plot(idx_stats)
