import django
from django.conf import settings
from django.template.loader import get_template
import sys, os

RESDIR = "./test/res"


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


def collect_res():

    resdict={}
    for file in os.listdir(RESDIR):
        if file.endswith(".html"):
            fname = os.path.splitext(file)[0]
            fpath =os.path.join(RESDIR, file)
            resdict[fname] = fpath

    return resdict


if __name__  == "__main__":

    resdict =collect_res()
    template = sys.argv[1]

    html = render(resdict, template)
    print(html)


