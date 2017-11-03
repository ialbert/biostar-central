from biostar.engine.settings import *
from .recipes_secrets import *

DEBUG = False

WSGI_APPLICATION = 'conf.recipes.recipes_wsgi.application'

ALLOWED_HOSTS = ['bioinformatics.recipes', 'www.bioinformatics.recipes']

SITE_HEADER =  "Bioinformatics Recipes"