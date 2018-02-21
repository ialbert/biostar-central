import hjson
import os

import django
from django.conf import settings
from django.template import Template, Context
from django.template import loader

__DIR = os.path.dirname(__file__)
__TEMPLATES = os.path.join(__DIR, "templates")


def setup():
    '''
    Specify template engine and template directory.
    Pass the settings variable to configure.
    '''

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': [__TEMPLATES],
            'OPTIONS': {
                'string_if_invalid': "** MISSING **",
                'libraries': {
                    'plotify': 'biostar.tools.plotify',
                    'igv': 'biostar.tools.igv.igv_tags',
                    'iobio': 'biostar.tools.iobio.iobio_tags',
                },
            },
        }
    ]
    settings.configure(TEMPLATES=TEMPLATES)
    django.setup()
    return


def read_data(fname):
    return hjson.load(open(fname))


def read_template(fname):
    return open(fname).read()


def render_file(data, fname):
    template_txt = open(fname).read()
    return render_string(data, template_txt)


def render_string(data, template_txt):
    spec = hjson.loads(data)
    return render_data(spec, template_txt)


def render_data(context, template_txt):
    setup()
    template = Template(template_txt)
    ctx = Context(context)
    html = template.render(ctx)
    return html


def render_template(data, name):
    '''
    Attempts to find and load a template by name
    via the django template loading mechanism.
    '''
    setup()
    template = loader.get_template(name)
    result = template.render(data)
    return result


if __name__ == "__main__":

    path1 = '/Users/ialbert/edu/24/bwa1.bam'
    name1 = os.path.basename(path1)

    path2 = '/Users/ialbert/edu/24/bwa2.bam'
    name2 = os.path.basename(path1)


    bams = [ (path1, name1),
             (path2, name2)
    ]

    data = dict(bams=bams)

    name = "igv/igv.xml"

    html = render_template(data, name)

    print(html)
