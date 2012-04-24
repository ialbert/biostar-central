#
# import from the main settings then override some of them
#
from main.settings import *

# set the domain administrators
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

# set the site url
SITE_DOMAIN = 'localhost:8080'

# set the secret key for the site
SECRET_KEY = 'my-secret key goes here'

# set your google tracker
GOOGLE_TRACKER = ""

# database setup
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': 'test-database',                  
        'USER': 'someuser',                      
        'PASSWORD': 'somepassword',                  
        'HOST': 'somehost',                      
        'PORT': 'somepost',                      
    }
}

# this sets wether to allow test logins via selenium
ALLOW_SELENIUM_TEST_LOGIN = False

# set up the email provider for your site
EMAIL_HOST = 'smtp.domain.org'
EMAIL_HOST_USER = 'user'
EMAIL_HOST_PASSWORD = 'password'
DEFAULT_FROM_EMAIL = 'email@domain.org'
SERVER_EMAIL = 'email@domain.org'