from biostar.engine.settings import *
from .recipes_secrets import *

DEBUG = False

WSGI_APPLICATION = 'conf.main.main_wsgi.application'

ALLOWED_HOSTS = ['bioinformatics.recipes', 'www.bioinformatics.recipes']
