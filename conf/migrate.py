#
# import from the main settings then override some of the values
#
from main.settings import *

# database setup
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': 'biostar-export',
        'USER': 'ialbert',
        'PASSWORD': '',
        'HOST': '',                      
        'PORT': '',                      
    }
}

