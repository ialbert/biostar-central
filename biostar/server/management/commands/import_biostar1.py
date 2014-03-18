from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from django.utils.timezone import utc
from django.conf import settings
from django.contrib.auth import get_user_model
import os, csv, datetime
from django.utils.dateparse import parse_datetime
from django.utils import timezone, encoding

from django.db import transaction
from itertools import *
from django.db.models import signals


def path_join(*args):
    return os.path.abspath(os.path.join(*args))

# Obtain the user model
User = get_user_model()

BATCH_SIZE = 100

USER_TYPE_MAP = {
    "New": User.NEW_USER, "Member": User.USER, "Blog": User.BLOG,
    "Moderator": User.MODERATOR, "Administrator": User.ADMIN
}

USER_STATUS_MAP = {
    "Active": User.TRUSTED, "Suspended": User.SUSPENDED, "Banned": User.BANNED,
}



MIGRATE_DIR = os.environ["BIOSTAR_MIGRATE_DIR"]
MIGRATE_DIR = os.path.expanduser(MIGRATE_DIR)

def get(data, attr, func=encoding.smart_unicode):
    value = data.get(attr, '').strip()
    try:
        value = func(value)
    except Exception, exc:
        raise Exception(value)
    return value

def localize_time(text):
    naive = parse_datetime(text)
    local = timezone.make_aware(naive, timezone=utc)
    return local

def get_post(row, users, klass):

    POST_TYPE_MAP = {
        "Question": klass.QUESTION, "Answer": klass.ANSWER,
        "Comment": klass.COMMENT, "Job": klass.JOB, "Blog": klass.BLOG,
    }


    uid = get(row, 'id', func=int)
    root_id = get(row, 'root_id', func=int)
    parent_id = get(row, 'parent_id', func=int)

    title = get(row, 'title').title()
    title = title[:140]
    tag_val = get(row, 'tag_val').strip()

    author_id = get(row, 'author_id', func=int)
    author = users.get(author_id)

    if not author:
        print("*** author %s not found for post %s" % (author_id, uid))
        return None

    post_type = get(row, 'post_type')
    post_type = POST_TYPE_MAP.get(post_type, klass.FORUM)
    post_status = klass.OPEN if get(row, 'post_status') == "Open" else klass.CLOSED

    post = klass(id=uid, title=title, author=author, lastedit_user=author,
                parent_id=parent_id, root_id=root_id)

    post.status = post_status
    post.type = post_type
    post.tag_val = ",".join(tag_val.split())
    post.creation_date = localize_time(get(row, 'creation_date'))
    post.lastedit_date = localize_time(get(row, 'lastedit_date'))
    post.view_count = get(row, "views", func=int)
    post.reply_count = get(row, "answer_count", func=int)
    post.book_count = get(row, "book_count", func=int)
    post.thread_score = get(row, "full_score", func=int)
    post.vote_count = get(row, "score", func=int)

    post_file = path_join(MIGRATE_DIR, 'posts', str(uid))
    post.content = file(post_file, 'rt').read()

    return post

def to_unicode(obj, encoding='utf-8'):
    if isinstance(obj, basestring):
        if not isinstance(obj, unicode):
            obj = unicode(obj, encoding)
    return obj


class Command(BaseCommand):
    help = 'migrate data from Biostar 1.*'

    option_list = BaseCommand.option_list + (
        make_option("-u", '--users', action="store_true", dest='users', default=False, help='import users'),
        make_option("-p", '--posts', action="store_true", dest='posts', default=False, help='import posts'),
        make_option("-x", '--votes', action="store_true", dest='votes', default=False, help='import votes'),
        make_option("-t", '--tags', action="store_true", dest='tags', default=False, help='auto tag'),
    )

    def auto_tag(self, fname):
        pass


    def migrate_posts(self, fname):
        from biostar.apps.posts.models import Post, Subscription
        from biostar.apps.messages.models import Message

        log = self.stdout.write

        # Disconnect the message creation signal,
        # it would generate way too many messages
        signals.post_save.disconnect(dispatch_uid="post-create-messages")

        Post.objects.all().delete()

        users = dict((u.id, u) for u in User.objects.all())

        log("migrating posts from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        for i, row in enumerate(stream):
            title = to_unicode(row['title'])

            log("migrating %s: %s" % (i, title))
            post = get_post(row, users, klass=Post)

            if not post:
                continue
            post.save()

            # TODO migrate only tags with high counts
            post.add_tags(post.tag_val)

            #log("migrated %s" % post)

        log("migrated %s posts" % Post.objects.all().count())
        log("created %s subscriptions" % Subscription.objects.all().count())
        log("created %s messages" % Message.objects.all().count())

    def migrate_users(self, fname):

        #User.objects.all().delete()
        log = self.stdout.write

        log("migrating users from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        seen = set()

        for row in stream:
            uid = int(get(row, 'id'))

            # Skip the first user. It is the default admin.
            if uid == 1:
                continue

            email = get(row, 'email')
            name = get(row, 'display_name')
            score = get(row, 'score')
            scholar = get(row, 'scholar')
            location = get(row, 'location')
            website = get(row, 'website')
            is_active = get(row, 'status') == "Active"
            user_type = USER_TYPE_MAP[get(row, 'type')]
            date_joined = get(row, 'date_joined')
            last_visited = get(row, 'last_visited')

            # Start populating the user.
            user = User(id=uid, email=email)
            user.email = email or "%s@xyz.xyz" % uid
            user.name = name
            user.score = score
            user.type = user_type

            # Original email were not required to be unique.
            if user.email in seen:
                user.email = "%s@biostars.org" % uid

            # Adds the email to the seen bucket.
            seen.add(user.email)

            user.is_active = is_active
            user.save()

            # Populate the profile
            prof = user.profile
            prof.website = website
            prof.scholar = scholar
            prof.location = location
            prof.date_joined = localize_time(date_joined)
            prof.last_login = localize_time(last_visited)
            about_me_file = path_join(MIGRATE_DIR, 'about_me', str(uid))
            prof.about_me = file(about_me_file, 'rt').read()
            prof.save()

            log("migrated user %s:%s" % (user.id, user.email))

        if settings.DEBUG:
            for id in (2, 10,):
                try:
                    # We use this during debugging to make it easy to log in as someone else
                    bot = User.objects.get(id=id)
                    log(
                        "updated user %s with email=%s, name=%s, password=SECRET_KEY," % (bot.id, bot.email, bot.name))
                    bot.set_password(settings.SECRET_KEY)
                    bot.save()
                except Exception, exc:
                    pass

        log("migrated %s users" % User.objects.all().count())


    def migrate_votes(self, fname):
        log = self.stdout.write

        from biostar.apps.posts.models import Vote

        VOTE_TYPE_MAP = {
            "Upvote": Vote.UP, "Downvote": Vote.DOWN,
            "Accept": Vote.ACCEPT, "Bookmark": Vote.BOOKMARK,
        }


        Vote.objects.all().delete()

        log("migrating votes from %s" % fname)

        user_map = dict((u.id, u) for u in User.objects.all())

        stream = csv.DictReader(file(fname), delimiter=b'\t')

        def vote_generator():
            for i, row in enumerate(stream):
                post_id = int(get(row, 'post_id'))
                author_id = int(get(row, 'author_id'))
                author = user_map.get(author_id)

                if not author:
                    log("*** author %s not found" % author_id)
                    continue

                vote_type = get(row, 'vote_type')
                vote_type = VOTE_TYPE_MAP[vote_type]
                vote_date = get(row, 'vote_date')
                vote_date = localize_time(vote_date)

                # Create the vote.
                vote = Vote(author=author, post_id=post_id, type=vote_type, date=vote_date)
                log("*** adding votes %s" % i)
                yield vote

        # Insert votes in batch. Bypasses the signals!
        Vote.objects.bulk_create(vote_generator(), batch_size=500)

        log("migrated %s votes" % Vote.objects.all().count())

    def handle(self, *args, **options):

        if options['users']:
            fname = path_join(MIGRATE_DIR, "users.txt")
            self.migrate_users(fname)

        if options['posts']:
            fname = path_join(MIGRATE_DIR, "posts.txt")
            self.migrate_posts(fname)

        if options['votes']:
            fname = path_join(MIGRATE_DIR, "votes.txt")
            self.migrate_votes(fname)

        if options['tags']:
            fname = path_join(MIGRATE_DIR, "posts.txt")
            self.auto_tag(fname)
