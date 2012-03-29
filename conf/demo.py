"""
How to change the settings for a different site
"""

# import from the main settings then override some of them
from main.settings import *

# set the domain administrators
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

# set the site url
SITE_DOMAIN = 'localhost:8080'

# set the secret key for the site
SECRET_KEY = 'my-secret key goes here'

# your google tracker
GOOGLE_TRACKER = "UA-300333-13"

# a local postgresql database
DATABASE_NAME = 'biostar-test-database'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': DATABASE_NAME,                  
        'USER': '',                      
        'PASSWORD': '',                  
        'HOST': '',                      
        'PORT': '',                      
    }
}

# set up the email provider for your site
EMAIL_HOST = 'smtp.domain.org'
EMAIL_HOST_USER = 'user'
EMAIL_HOST_PASSWORD = 'password'
DEFAULT_FROM_EMAIL = 'email@domain.org'
SERVER_EMAIL = 'email@domain.org'