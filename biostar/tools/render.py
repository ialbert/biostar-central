import django
from django.conf import settings
from django.template import Template,Context
import os, hjson


def set_env():
    '''
    Specify template engine and template directory.
    Pass the settings variable to configure.
    '''

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': '',
            'OPTIONS': {
                'string_if_invalid': "***MISSING**"
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


def render_data(spec, template_txt):
    set_env()
    template = Template(template_txt)
    context = Context(spec)
    html = template.render(context)
    return html


if __name__ == "__main__":

    config = "./templates/metabarcode_qc/metabarcode_spec.hjson"
    template = "./templates/metabarcode_qc/metabarcode_makefile.html"

    html = render_file(config, template)
    print(html)





