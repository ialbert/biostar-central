from biostar.tools.plotify import ChartParams,BarchartParams
from biostar.tools.render import render_template
import sys, csv , re
import csv


def create_factors(fname, selected_reads):
    '''
    factor is calculated as total_reads/selected_reads.
    returns a dictionary with factor calculated for each sample.
    '''
    stream = open(fname).readlines()

    factors = dict()

    for data in stream:
        if data.strip().isdigit():
            count = int(data.strip())
            factors[sample] = float(count) / float(selected_reads)
        sample = data.strip()

    return factors


def get_factor(factors, sample):
    if sample in factors:
        return factors[sample]
    return 1


def parse_idxstats(fname, factors):
    '''
    parse samtools idxstats and multiply normalize the mapped counts as
    normalized = (mapped_counts/scaffold_length)*factor
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

        if name == "*" :
            continue

        # normalize by scaffold length
        len_norm = float(mapped)/float(size)
        mapped_cov = len_norm * get_factor(factors, sample)
        elem.setdefault(name, []).append((int(size), float(mapped_cov)))

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

    return rows


def plot(data1):

    # Plot 2
    p1 = BarchartParams()
    p1.type = 'BarChart'
    p1.header = data1[0]
    p1.data = data1[1:]
    p1.options = '''    
            title: 'Mapped counts per chromosome.',
            chartArea: {  width: "50%" }
        '''

    # This is the context.
    data = dict(p1=p1)

    name = "scaffold.html"

    html = render_template(data, name)
    print(html)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create barcharts from samtools idxstats output.')

    parser.add_argument('--mapped', help='File with idxstats output ')
    parser.add_argument('--total', help='File with total read counts ')
    parser.add_argument('--selected', help='No. of reads subselected for mapping')

    args = parser.parse_args()

    # Read the arguments
    mapped = args.mapped
    total = args.total
    selected = args.selected

    # Normalization
    #(mapped_counts/scaffold_length)*(total_read_count/selected_read_count)

    # Calculate factors to normalize with
    factors = create_factors(total, selected)

    # Parse idxstats into a dictionary after normalizing mapped counts.
    store = parse_idxstats(mapped, factors)

    # Format the idxstats dictionary into a list of lists.
    idx_stats = format_idxstats(store)

    # Plot barchart.
    plot(idx_stats)

