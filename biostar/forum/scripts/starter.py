# Starter script template for writing programs that
# interact with the django site
import os
import sys

# Biostar home directory
PATH = "~/app/biostar-central"

PATH = os.path.expanduser(PATH)

sys.path.append(PATH)

import django

from django.conf import settings

#os.environ.setdefault("DJANGO_SETTINGS_MODULE", "biostar.forum.settings")

# Import the settings file
from biostar.forum import settings as app_settings

#print (os.environ.get('DJANGO_SETTINGS_MODULE', "*** NOT FOUND ***") )

# Apply settings onto the django settings module.
settings.configure(default_settings=app_settings, DEBUG=True, LOGGING_CONFIG=None, FORCE_SCRIPT_NAME=None)

# Bootstrap the Django framework.
django.setup()

# You can now import and user models in the app.
from biostar.forum.models import User



def main():
    users = User.objects.all().count()
    print(users)
    pass


if __name__ == '__main__':
    main()
