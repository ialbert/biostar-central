#
# import from the main settings then override some of the values
#
from main.settings import *

# set the domain administrators
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

# set the site url
SITE_DOMAIN = 'localhost:8080'

# set the cookie domain as needed
SESSION_COOKIE_DOMAIN = ".biostars.org"

# set the secret key for the site
SECRET_KEY = 'my-secret key goes here'

# set your google tracker
GOOGLE_TRACKER = ""

# this is the authentication module and function in it
EXTERNAL_AUTHENTICATOR_FUNC = 'dummy_auth'

MIDDLEWARE_CLASSES.extend([
    'main.extauth.ExternalAuthenticator',
])

# database setup
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': 'test-database',                  
        'USER': 'someuser',                      
        'PASSWORD': 'somepassword',                  
        'HOST': '',                      
        'PORT': '',                      
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

patt = """Galaxy username <b>%(display_name)s</b>\n
Tool name: <b>%(tool_name)s</b>
Tool version:  <b>%(tool_version)s</b>
Tool id: <b>%(tool_id)s</b>
Tags: <b>%(tags)s</b>"""

EXTERNAL_AUTHENICATION = {
    "TEST-KEY" : ("abcd", patt),
}

