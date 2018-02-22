import sys,csv
from biostar.tools.plotify import ChartParams
from biostar.tools.render import render_template


def parse_classification(fname):
    """parses centrifuge classification report file."""

    items = []
    # Read centrifuge report file.
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            species = row['name']
            reads = row['numUniqueReads']
            items.append((species, reads))

    # Sort by counts.
    # items.sort(key=lambda x: int(x[1]))

    # Sort by species name.
    items.sort(key=lambda x: x[0])

    return items


def plot(data):
    # Plot 1

    p1 = ChartParams()
    p1.type = 'BarChart'
    p1.data = data
    p1.xlabel = "Species"
    p1.ylabel = "Read counts"
    p1.options = '''    
            title: 'Unique reads in species',
            legend: {position: 'none'},
        '''

    # Plot 2

    p2 = ChartParams()
    p2.type = 'PieChart'
    p2.data = data
    p2.xlabel = "Species"
    p2.ylabel = "Read counts"
    p2.options = '''
            title: 'Unique reads in species',
            chartArea:{left:20,top:0,width:'80%',height:'75%'},
            legend:{position:'left'}

        '''

    # Plot3

    locations = [(41.015908, -77.53246, 'Lamar PA')]
    p3 = ChartParams()
    p3.type = 'Map'
    p3.data = locations
    p3.options = '''
            zoomLevel: 6,
            showTooltip: true,
            showInfoWindow: true

        '''

    # This is the context.
    data = dict(p1=p1, p2=p2, p3=p3)

    name = "classify.html"

    html = render_template(data, name)

    print(html)
    #with open('index.html', 'wt') as fp:
    #    fp.write(html)


if __name__ == '__main__':

    fname = sys.argv[1]
    #fname = "report.txt"
    results = parse_classification(fname)
    plot(results)



