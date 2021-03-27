#
# Starter script template for writing programs that interact with the django site
#

import os, sys, plac
from itertools import *
import django
from django.conf import settings
from datetime import timedelta

# Biostar home directory
PATH = "~/biostar-central"
PATH = os.path.expanduser(PATH)

# Add biostar-central to the Python path.
sys.path.append(PATH)

MODULE = "DJANGO_SETTINGS_MODULE"

# If module not set, use a sane default.
if MODULE not in os.environ:
    os.environ.setdefault(MODULE, "biostar.forum.settings")

# A reminder of what the settings are
#print(f"# {MODULE}={os.environ[MODULE]}")

# Bootstrap the Django framework.
django.setup()

from biostar.forum import util, auth
from biostar.forum.models import User, Profile, Post, Log, Vote


def get_admin():
    admin = User.objects.filter(is_superuser=True).order_by("pk").first()
    return admin

def get_posts(user):
    posts = Post.objects.filter(author=user)
    return posts

def recently(days=2):
    return util.now() - timedelta(days=days)


@plac.opt("limit", "limit ", abbrev="L")
@plac.opt("days", "deletes the users", abbrev="d")
@plac.flg("delete", "deletes the users", abbrev='D')
def main(delete=False, days=2, limit=10):
    admin = get_admin()
    
    posts = Post.objects.filter(spam=Post.SPAM, author__profile__state=Profile.NEW).order_by("-pk")
    posts = islice(posts, limit)

    users = {}
    for post in posts:
        #print (post.title)
        users[post.author.id] = post.author
    
    u_count = s_count = 0

    # Check recent users for spamming.
    for user in users.values():
        
        # How many posts did the user make.
        posts = get_posts(user)
        post_count = posts.exclude(spam=Post.SPAM).count()

        # How many spam posts did the user make.
        spam_count = posts.filter(spam=Post.SPAM).count()

        if spam_count > post_count:
            u_count += 1
            s_count += spam_count
            text = f"user={user.profile.name} spam_count={spam_count}"
            #print (text)
            if delete:
                user.delete()
                #print("deleted")

    if delete and s_count:
        msg = f"spam cleanup, removed {u_count} spammers and {s_count} spam posts"
        auth.db_logger(user=admin, text=msg)

if __name__ == "__main__":
    plac.call(main)
