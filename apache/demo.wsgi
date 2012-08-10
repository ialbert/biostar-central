"""
This setup is specific to our hosting company
"""
import os, sys

# customize this
DJANGO_SETTINGS_MODULE = 'conf.demo'

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

# the directory that this file is located in
__CURR_DIR = path(os.path.dirname(__file__))

# the root folder that contains dependencies
__ROOT = path('..', __CURR_DIR)

extras = [
        path(__ROOT),
        path(__ROOT, 'libs'),
        path(__ROOT, 'libs', 'libraries.zip'),
]

# add the libraries to the import path
sys.path.extend( extras )

from django.core.handlers.wsgi import WSGIHandler
os.environ['DJANGO_SETTINGS_MODULE'] = DJANGO_SETTINGS_MODULE
application = WSGIHandler()
