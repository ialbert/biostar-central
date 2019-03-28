
from biostar.settings import *

USERS_DATABASE_NAME = join(BASE_DIR, '..', 'export', 'database', 'users.db')


DATABASES.update(
    {'users':
        {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': USERS_DATABASE_NAME,
        }
    } )
