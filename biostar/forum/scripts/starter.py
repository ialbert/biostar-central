# Starter script template for writing programs that
# interact with the django site
import sys, os

# Biostar home directory
DIR="/home/ialbert/app/biostar-central"

sys.path.append(DIR)

import django
from django.conf import settings


# Import the settings file
from biostar.forum import settings as app_settings

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "biostar.forum.settings")

# Apply settings onto the django settings module.
#settings.configure(default_settings=app_settings, DEBUG=True, LOGGING_CONFIG=None)


# Bootstrap the Django framework.
django.setup()

# You can now import and user models in the app.
from biostar.forum.models import User, Post

def main():
    users = User.objects.filter().all()
    print(users)
    pass


if __name__ == '__main__':
    main()