from biostar.tools.plotify import ChartParams, TableParams
from biostar.tools.render import render_template
import sys, csv , re


class Table():
    def __init__(self):
        self.rows =[]
        self.colums = []


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
        stats.append((tag.strip(), value.strip()))

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


def clean_accession(chrom):

    patt = r'gi\|[0-9]{1,}\|[a-z]{1,}\|'
    chrom = re.sub(patt, '', chrom)
    chrom = chrom.replace('|', '')
    return chrom


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def parse_table(fname):
    '''
    convert chromosome classification table to goolge table format.
    '''

    def get_types(first_line):
        type_list = []
        for cols in first_line:
            if is_int(cols):
                type_list.append('number')
            else:
                type_list.append('string')
        return type_list

    stream = open(fname).readlines()

    # Get columns.
    col_names = stream[0].split("\t")
    col_names[-1] = col_names[-1].strip()

    first_line = stream[1].split("\t")
    col_types = get_types(first_line)

    columns = list(zip(col_types,col_names))

    # Get rows.
    rows = []
    data = stream[1:]
    for line in data:
        arr = line.split("\t")
        row=[]
        for item in arr:
            item = item.strip()
            if is_int(item):
                item = int(item)

            else:
                item = clean_accession(item)
                item = "'{0}'".format(item)
            row.append(item)
        #print(row)
        rows.append(row)
    return columns, rows


def plot(data1, data2 , data3):

    # Plot 1
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
            title: '',
            legend: {position: 'none'},
            colors : ['green']
        '''

    # Plot 2
    p2 = ChartParams()
    p2.type = 'BarChart'
    p2.data = data2
    p2.xlabel = "Chromosomes"
    p2.ylabel = "Read counts"
    p2.options = '''    
            title: '',
            legend: {position: 'none'},
            colors: ['green']
        '''

    # Plot table.
    t1 = TableParams()
    t1.rows = data3.rows
    t1.columns = data3.columns

    t1.options = '''    
               width : '100%',
               height : '50%'
           '''

    # This is the context.
    data = dict(p1=p1, p2=p2, t1=t1)

    name = "bwa_classify.html"

    html = render_template(data, name)
    print(html)

    #with open("index.html", "w") as fh:
    #    fh. write(html)


if __name__ == '__main__':

    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    fname3 = sys.argv[3]

    #fname1 = "alignment_stats.txt"
    #fname2 = "chrom_mapping.txt"
    #fname3 = "classification_table.txt"

    flag_stats = parse_flagstats(fname1)
    idx_stats = parse_idxstats(fname2)
    columns, rows = parse_table(fname3)
    table = Table()
    table.rows = rows
    table.columns = columns

    plot(flag_stats, idx_stats, table)
