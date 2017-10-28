from biostar.tools.render import render_template
from django import template
from itertools import count
register = template.Library()

@register.inclusion_tag('charts/chart_init.html')
def chart_init():
    context = dict()
    return context

@register.inclusion_tag('charts/chart_maker.js')
def chart_maker(params):
    context = dict(params=params)
    return context

@register.inclusion_tag('charts/table_maker.js')
def table_maker(params):
    context = dict(params=params)
    return context


counter = count(100)


class ChartParams():
    def __init__(self):
        self.data = []
        self.type = 'Histogram'
        self.xlabel = self.ylabel = 'label'
        self.options = ''
        self.count_id = next(counter)
        self.chart_id = f'chart_{self.count_id}'


class TableParams():
    def __init__(self):
        self.type = 'Table'
        self.rows = []
        self.columns = []
        self.options = ''
        self.count_id = next(counter)
        self.table_id = f'table_{self.count_id}'


def demo():
    from random import randint

    def make(x):
        value = randint(1, 5)
        return f'A{value}', value


    p1 = ChartParams()
    p1.data = map(make, range(50))

    p1.xlabel = "Size 1"
    p1.ylabel = "Value 1"

    p1.options = '''    
        title: 'Histogram Plot 1',
        legend: {position: 'none'},
    '''

    p2 = ChartParams()
    p2.type = 'BarChart'
    p2.data = map(make, range(15))
    p2.xlabel = "Size 2"
    p2.ylabel = "Value 2"

    p2.options = '''    
           title: 'Histogram Plot 2' ,
           legend: {position: 'none'},
       '''

    # Plot table.

    t1 = TableParams()
    rows = [
        ["'Mike'", 10000, 'true'],
        ["'Jim'", 8000, 'false'],
        ["'Alice'", 12500, 'true'],
        ["'Bob'", 7000, 'false']
    ]
    columns = [
        ('string', 'Name'),
        ('number', 'Salary'),
        ('boolean', 'Full Time Employee')
    ]
    t1.columns = columns
    t1.rows = rows
    t1.options = '''    
            title: 'Table 1',
            width: '50%',
            height: '50%'
        '''

    # This is the context.
    data = dict(p1=p1, p2=p2, t1=t1)

    name = "chart_demo.html"

    html = render_template(data, name)

    with open('index.html', 'wt') as fp:
        fp.write(html)


if __name__ == '__main__':
    demo()
