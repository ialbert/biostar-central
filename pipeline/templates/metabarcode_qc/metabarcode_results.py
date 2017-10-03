#from ...make import render
import django
from django.conf import settings
from django.template.loader import get_template
import sys, os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_DIR = CURR_DIR


def set_env():
    '''
    Specify template engine and template directory.
    Pass the settings variable to configure.
    '''

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': ['./'],
        }
    ]
    settings.configure(TEMPLATES=TEMPLATES)
    return


def render(data,template):
    set_env()
    django.setup()
    temp = get_template(template)
    html = temp.render(data)
    return html


def collect_res(resdir):

    resdict={}
    for file in os.listdir(resdir):
        if file.endswith(".html"):
            fname = os.path.splitext(file)[0]
            fpath =os.path.join(resdir, file)
            resdict[fname] = fpath

    return resdict


if __name__  == "__main__":

    resdir = sys.argv[1]
    template = sys.argv[2]
    results =collect_res(resdir)
    context={'results' : results}
    html = render(context, template)


