from biostar.engine.settings import *
from .main_secrets import *

DEBUG = False

WSGI_APPLICATION = 'conf.main.main_wsgi.application'

ALLOWED_HOSTS = ['metabarcode.com', 'www.metabarcode.com']
