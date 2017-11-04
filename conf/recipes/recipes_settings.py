from biostar.engine.settings import *
try:
    from .recipes_secrets import *
except ImportError as exc:
    print("No recipes.recipes_secrets module could be found")
DEBUG = False

WSGI_APPLICATION = 'conf.recipes.recipes_wsgi.application'

ALLOWED_HOSTS = ['bioinformatics.recipes', 'www.bioinformatics.recipes']

SITE_HEADER =  "Bioinformatics Recipes"