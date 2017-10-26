import django
from django.conf import settings
from django.template import Template,Context
from django.template import loader

import os, hjson
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
            'DIRS': [ __TEMPLATES ],
            'OPTIONS': {
                'string_if_invalid': "** MISSING **",
                'libraries': {
                    'plotify': 'biostar.tools.plotify',
                },
            },
        }
    ]
    settings.configure(TEMPLATES=TEMPLATES)
    django.setup()
    return


def render_file(data, template_file):
    template_txt = open(template_file).read()
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
    Attempts to find and load a template by name from
    via the django template loading mechanism.
    '''
    setup()
    template = loader.get_template(name)
    result = template.render(data)
    return result

if __name__ == "__main__":

    data = dict(name="World")
    name = "hello.html"

    html = render_template(data, name)

    print(html)





