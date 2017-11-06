from biostar.settings import *

DEBUG = True

WSGI_APPLICATION = 'conf.devel.devel_wsgi.application'

ALLOWED_HOSTS = ['localhost', 'www.lvh.me']

try:
    from .devel_secrets import *
except ImportError as exc:
    print("No devel_secrets module could be imported")
