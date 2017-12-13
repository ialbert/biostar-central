from biostar.settings import *

DEBUG = True

WSGI_APPLICATION = 'conf.email.email_wsgi.application'

SITE_ID = 1
SITE_DOMAIN = "psu.bioinformatics.recipes"
SITE_NAME = "PSU Bioinformatics Recipes"

HTTP_PORT = ''
PROTOCOL = 'https'

ALLOWED_HOSTS = [SITE_DOMAIN]

SITE_HEADER = '<img class="ui middle aligned small image" src="/static/images/shield.png"> <span>Penn State Bioinformatics Recipes</span>'

try:
    from .email_secrets import *
except ImportError as exc:
    print("*** No email_secrets module could be imported")
