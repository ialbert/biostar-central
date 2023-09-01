# Psycopg based database
from biostar.forum.settings import *

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": "biostardb",
        "USER": "test",
        "PASSWORD": "test",
    }
}
