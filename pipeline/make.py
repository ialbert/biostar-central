import django
from django.conf import settings
from django.template.loader import get_template
import sys, os, json, hjson
from const import *


def path_join(*args):
    return os.path.abspath(os.path.join(*args))


CURR_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_DIR = CURR_DIR
#TEMPLATE_DIR = path_join(CURR_DIR, "templates")


def set_env():
    '''
    Specify template engine and template directory.
    Pass the settings variable to configure.
    '''

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': [TEMPLATE_DIR],
        }
    ]
    settings.configure(TEMPLATES=TEMPLATES)
    return


def render(data, template):
    set_env()
    django.setup()
    temp = get_template(template)
    html = temp.render(data)
    return html


def name_val_pair(element):
    '''
    input is a dictionary object and it returns a name, value pair.
    '''

    try:
        value = element['default'] if not element['value'] else element['value']
    except KeyError:
        print("{0} raises an error".format(element['name']) )
    return element['name'], value


def parse_configs(configs):
    '''
    parse config file and returns a dictionary to be rendered.
    '''

    specs = dict()
    for element in configs:
        if not element['name'] == "analysis_spec":
            name, value = name_val_pair(element)
            specs[name] = value
    return specs


def load_specs(fname):
    try:
        data = open(fname).read()
    except OSError as e:
        print("File not found: {0}".format(e))
        sys.stderr.close()
    #data = json.loads(data)
    data = hjson.loads(data)
    return data


def run():
    config = sys.argv[1]
    template = sys.argv[2]

    #config = "./templates/metabarcode_qc/metabarcode_spec.hjson"
    #template = "./templates/metabarcode_qc/metabarcode_makefile.html"


    # loads sepc file.
    configs = load_specs(config)

    # parse configs to render in template.
    specs = parse_configs(configs)

    context = {'specs': specs, 'dirinfo': DIR_INFO, 'jobid': 'job0' }

    # render template.
    html = render(context, template)
    print(html)


if __name__ == "__main__":
    run()
