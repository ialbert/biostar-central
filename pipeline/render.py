import django
from django.conf import settings
from django.template.loader import get_template
import json
import sys,os


def set_env():
    '''
    Specify template engine and template directory.
    Pass the settings variable to configure.
    '''

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': ['./templates'],
        }
    ]
    settings.configure(TEMPLATES=TEMPLATES)


def dict2json(datadict):
    return json.dumps(datadict)


def metadata_loader(fname):
    '''
    reads json file/string and returns the data dictionary
    :param fname:
    :return: data dictionary.
    '''
    if os.path.isfile(fname):
        try:
            data = open(fname).read()

        except OSError as e:
            print("File not found: {0}".format(e))
            sys.stderr.close()
    else:
        data = fname
    data = json.loads(data)
    return data


def render_template(data, template_file):
    set_env()
    django.setup()
    data = metadata_loader(data)
    temp = get_template(template_file)
    html = temp.render(data)
    return html


if __name__ == '__main__':

    metadata = sys.argv[1]  # metadata as a json string or json file
    template = sys.argv[2]  # template file

    rendered = render_template(metadata, template)
    print(rendered)





