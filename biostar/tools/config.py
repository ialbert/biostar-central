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

try:

    # Bootstrap the Django framework.
    django.setup()

    # Get the settings.
    from django.conf import settings

    # Import Biostar specific modules.
    from biostar.forum import util, auth
    from biostar.forum.models import User, Profile, Post, Log, Vote
except Exception as exc:
    print(f"*** module: {__name__}")
    print(f"*** settings: {MODULE_NAME}={os.environ[MODULE_NAME]}")
    print(f"*** error: {exc}")
    sys.exit(1)

def past_date(days=2, minutes=0, seconds=0):
    return util.now() - timedelta(days=days, minutes=minutes, seconds=seconds)


def get_first_admin():
    """
    Returns the first admin among the users.
    """
    admin = User.objects.filter(is_superuser=True).order_by("pk").first()
    return admin


def get_posts(user):
    """
    Returns a query to posts by a user.
    """
    posts = Post.objects.filter(author=user)
    return posts
