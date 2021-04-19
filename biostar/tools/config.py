"""
Call this first in any tool. Bootstraps the django framework into a script.
"""
import os, sys
import django
from datetime import timedelta

# Module environment variable.
MODULE_NAME = "DJANGO_SETTINGS_MODULE"

# The default module to boot
DEFAULT_MODULE = "biostar.forum.settings"

# If module not set, use a sane default.
if MODULE_NAME not in os.environ:
    print(f"*** using the default settings module: {DEFAULT_MODULE}")
    os.environ.setdefault(MODULE_NAME, DEFAULT_MODULE)

# This is the module name we run with
module_name = os.environ[MODULE_NAME]
print(f"*** {MODULE_NAME}={module_name}")

try:
    # Get the settings.
    from django.conf import settings

    # Bootstrap the Django framework.
    django.setup()

    # Import Biostar specific modules.
    from biostar.forum import util, auth
    from biostar.forum.models import User, Profile, Post, Log, Vote
except Exception as exc:
    print(f"*** module={__name__}")
    print(f"*** {MODULE_NAME}={module_name}")
    print(f"*** error: {exc}")
    sys.exit(1)

def time_ago(days=2, minutes=0, seconds=0):
    return util.now() - timedelta(days=days, minutes=minutes, seconds=seconds)

def get_first_admin():
    """
    Returns the first admin among the users.
    """
    admin = User.objects.filter(is_superuser=True).order_by("pk").first()
    return admin


if __name__ == '__main__':
    print("*** running a configuration test")