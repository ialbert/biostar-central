__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os, re
from optparse import make_option
import random
import logging
from datetime import timedelta
from django.db.models import signals, Q
import string

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = "Runs quick patches over the data. Use it only if you know what you're doing."

    option_list = BaseCommand.option_list + (
        make_option('--users', dest='users', action='store_true', default=False, help='patches_users'),
        make_option('--bump', dest='bump', action='store_true', default=False, help='bumps a random post'),
        make_option('--bump_id', dest='bump_id', type=int, help='bumps a specific post'),
        make_option('--stuff', dest='stuff', action='store_true', default=False, help='runs stuff ...'),
        make_option('--tag', dest='tag', default="", help='tags post by matching a regex.Format regex:name'),
        make_option('--dry', dest='dry', action='store_true', default=False, help='dry run, sometimes applies ;-)'),
        make_option('--merge_users', dest='merge', metavar="FILE", default=False, help='merges users listed in a file, on per row: master alias1 alias2 ...'),
    )

    def handle(self, *args, **options):

        tag = options['tag']
        dry = options['dry']
        merge = options['merge']

        if tag:
            tagger(tag, dry)

        if options['stuff']:
            stuff()

        if options['users']:
            patch_users()

        if options['bump']:
            bump()

        if merge:
            merge_users(merge)

        pk = options['bump_id']
        if pk:
            bump(pk)



def post_patch():
    "One off tasks go here that just need a quick access to the data"
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User, Profile

    for post in Post.objects.all():
        post.html = post.content
        post.save()
    for prof in Profile.objects.all():
        prof.location = prof.location.strip()
        prof.save()


def stuff():
    "One off tasks go here that just need a quick access to the data"
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User, Profile
    from biostar.const import ALL_MESSAGES

    cond = Q(profile__message_prefs=ALL_MESSAGES)
    cond = Q(profile__tags__name__in=["galaxy"])
    users = User.objects.filter(cond)

    for user in users:
        print user.id, user.email


def tagger(pattern, dry):

    "One off tasks go here that just need a quick access to the data"
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User

    posts = Post.objects.filter(type__in=Post.TOP_LEVEL)
    patt, name = pattern.split(":")

    logger.info('%s -> %s' % (patt, name))

    patt = re.compile(patt, re.MULTILINE | re.IGNORECASE| re.DOTALL)
    for post in posts:
        try:
            hits = patt.search(post.content)
            if hits:
                logger.info(post.title)
                if not dry:
                    tag_val = "%s, %s" % (post.tag_val, name)
                    post.tag_val = tag_val
                    post.save()
                    post.add_tags(tag_val)
        except Exception, exc:
            logger.error("exception:'%s' while tagging %s: %s" % (exc, post.id, post.title))

def merge_users(fname):
    from biostar.apps.posts.models import Post, Subscription, ReplyToken
    from biostar.apps.posts.models import Vote
    from biostar.apps.users.models import User
    from biostar.apps.messages.models import Message, MessageBody
    from allauth.socialaccount.models import SocialAccount
    from allauth.account.models import EmailAddress

    print ("Merging users from %s" % fname)
    def create_pairs(line):
        elems = line.split()
        elems = map(string.strip, elems)
        first, rest = elems[0], elems[1:]
        if first in rest:
            msg = "master email among aliases: %s" % elems
            raise Exception(msg)
        return first, rest

    stream = map(string.strip, open(fname))
    stream = filter(None, stream)
    pairs = dict(map(create_pairs, stream))

    for key, values in pairs.items():

        master = User.objects.get(email=key)
        aliases = User.objects.filter(email__in=values).order_by('profile__date_joined')

        if not aliases:
            print("*** no matching aliases for master: %s, aliases: %s" % (key, ",".join(values)))
            continue

        # Keep the oldest dates for join and login
        master.profile.date_joined = aliases[0].profile.date_joined
        master.profile.last_login = aliases[0].profile.last_login
        master.profile.save()

        print("*** merging master: %s, aliases: %s" % (key, ",".join(values)))

        Post.objects.filter(author__in=aliases).update(author=master)
        Post.objects.filter(lastedit_user__in=aliases).update(lastedit_user=master)
        Vote.objects.filter(author__in=aliases).update(author=master)
        Message.objects.filter(user__in=aliases).update(user=master)
        MessageBody.objects.filter(author__in=aliases).update(author=master)

        # Migrate the social accounts
        SocialAccount.objects.filter(user__in=aliases).update(user=master)
        EmailAddress.objects.filter(user__in=aliases).update(user=master)

        # New score for the master
        score = sum( u.score for u in aliases ) + master.score
        User.objects.filter(pk=master.pk).update(score=score)

        # Update tokens
        ReplyToken.objects.filter(user__in=aliases).update(user=master)

        # Migrate subscriptions
        subs = Subscription.objects.filter(user__in=aliases)
        for sub in subs:
            # The master already has a subscription to the post.
            if Subscription.objects.filter(user=master, post=sub.post):
                sub.delete()
            else:
                sub.user = master
                sub.save()

        # Delete the users
        aliases.delete()


def patch_users():
    from biostar.apps.users.models import User, Profile
    from biostar.const import DEFAULT_MESSAGES
    #users = Profile.objects.all()
    #users.update(message_prefs=DEFAULT_MESSAGES)

def bump(pk=None):
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User

    from biostar.const import now

    if not pk:
        query = Post.objects.filter(type=Post.QUESTION, status=Post.OPEN)

        value = random.random()

        if value > 0.5:
            since = now() - timedelta(weeks=10)
            query = query.filter(reply_count=0, creation_date__gt=since)

        query = query.values_list("id")

        ids = [ p[0] for p in query ]

        pk = random.choice(ids)

    community = User.objects.get(pk=1)
    post = Post.objects.get(pk=pk)
    logger.info(post.title)

    if not post.is_toplevel:
        logger.warning("post is not at toplevel")

    post.lastedit_date = now()
    post.lastedit_user = community
    post.save()


