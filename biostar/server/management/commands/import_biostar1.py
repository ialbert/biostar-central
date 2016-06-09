from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from django.utils.timezone import utc, get_current_timezone
from django.conf import settings
from django.contrib.auth import get_user_model
import os, csv, datetime
from datetime import timedelta

from django.utils.dateparse import parse_datetime
from django.utils import timezone, encoding

from django.db import transaction
from itertools import *
from django.db.models import signals
from biostar.const import now

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

tz = get_current_timezone()

shift = timedelta(hours=6)

def get(data, attr, func=encoding.smart_unicode):
    value = data.get(attr, '').strip()
    try:
        value = func(value)
    except Exception as exc:
        raise Exception(value)
    return value

def localize_time(text):
    global shift
    naive = parse_datetime(text) + shift
    local = timezone.make_aware(naive, timezone=utc)
    return local

def get_post(row, users, klass):

    POST_TYPE_MAP = {
        "Question": klass.QUESTION, "Answer": klass.ANSWER,
        "Comment": klass.COMMENT, "Job": klass.JOB, "Blog": klass.BLOG,
        "Tool": klass.TOOL, "News": klass.NEWS,
        "Tutorial": klass.TUTORIAL,
    }

    POST_STATUS_MAP = {
        "Open": klass.OPEN,
        "Closed": klass.CLOSED,
        "Deleted": klass.DELETED,
    }

    uid = get(row, 'id', func=int)
    root_id = get(row, 'root_id', func=int)
    parent_id = get(row, 'parent_id', func=int)

    title = get(row, 'title').title()
    title = title[:200]
    tag_val = get(row, 'tag_val').strip()

    author_id = get(row, 'author_id', func=int)
    author = users.get(author_id)

    lastedit_user_id = get(row, 'lastedit_user', func=int)
    lastedit_user = users.get(lastedit_user_id)

    lastedit_user = lastedit_user or author

    if not author:
        print("*** author found for post %s" % (author_id, uid))
        return None

    post_type = get(row, 'post_type')

    post_type = POST_TYPE_MAP.get(post_type, klass.FORUM)

    if post_type == klass.TUTORIAL:
        tag_val += " tutorial"

    post_status = get(row, 'post_status')

    post_status = POST_STATUS_MAP[post_status]

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
    post.lastedit_user = lastedit_user

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
        make_option("-d", '--dir', dest='dir', default="~/tmp/biostar-migrate", help='import directory'),
    )

    def auto_tag(self, source, fname):
        pass


    def migrate_posts(self, source, fname):
        from biostar.server.models import disconnect_all
        from biostar.apps.posts.models import Post, Subscription
        from biostar.apps.messages.models import Message
        from biostar.apps.util import html

        log = self.stdout.write

        # Disconnect signals they will generate way too many messages
        disconnect_all()

        posts = [ p[0] for p in Post.objects.all().values_list("id") ]

        posts = set(posts)

        users = dict((u.id, u) for u in User.objects.all())

        log("migrating posts from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        for i, row in enumerate(stream):
            title = to_unicode(row['title'])
            uid = int(row['id'])
            url = row['url'].strip()

            # Skip existing posts
            if uid in posts:
                continue

            posts.add(uid)

            log("migrating post %s: %s" % (uid, title))
            post = get_post(row, users, klass=Post)

            if not post:
                log("skipped %s: %s" % (uid, title))
                continue

            # Read and add the post body.
            post_file = path_join(source, 'posts', str(post.id))
            post.content = file(post_file, 'rt').read()

            if url and post.type == Post.BLOG:
                # Will break out an not deal with Blogs in Biostar.
                continue
                # Link to external blog bosts.
                url_link = '<p><b><i class="fa fa-external-link-square"></i> Read full blogpost at <a href="%s">%s</a></b><p>' % (url, url[:45])
                url_link = to_unicode(url_link)
                content = to_unicode(post.content)
                post.content = url_link + content

            try:
                post.save()
            except Exception as exc:
                log('*** error inserting post %s' % post.id)
                log("*** %s" % exc)
                continue

            # TODO migrate only tags with high count
            post.add_tags(post.tag_val)


        log("migrated %s posts" % Post.objects.all().count())
        log("created %s subscriptions" % Subscription.objects.all().count())
        log("created %s messages" % Message.objects.all().count())

    def migrate_users(self, source, fname):

        #User.objects.all().delete()
        log = self.stdout.write

        log("migrating users from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        email_set, uid_seen = set(), set()

        users = dict((u.id, u) for u in User.objects.all())

        for row in stream:
            uid = int(get(row, 'id'))

            # The file may contain the same user multiple times
            # Caused by incremental dumping
            if uid in users or uid in uid_seen:
                continue

            uid_seen.add(uid)

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
            if user.email in email_set:
                user.email = "%s@biostars.org" % uid

            # Adds the email to the seen bucket.
            email_set.add(user.email)

            user.is_active = is_active
            user.save()

            # Populate the profile
            prof = user.profile
            prof.website = website
            prof.scholar = scholar
            prof.location = location
            prof.date_joined = localize_time(date_joined)
            prof.last_login = localize_time(last_visited)
            about_me_file = path_join(source, 'about_me', str(uid))
            prof.info = file(about_me_file, 'rt').read()
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
                except Exception as exc:
                    pass

        log("migrated %s users" % User.objects.all().count())


    def migrate_votes(self, source, fname):
        log = self.stdout.write

        from biostar.apps.posts.models import Post, Vote

        VOTE_TYPE_MAP = {
            "Upvote": Vote.UP, "Downvote": Vote.DOWN,
            "Accept": Vote.ACCEPT, "Bookmark": Vote.BOOKMARK,
        }

        Vote.objects.all().delete()

        posts = Post.objects.all().values_list('id')

        seen = set(p[0] for p in posts)

        log("loaded %s post ids" % len(seen) )

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

                if post_id not in seen:
                    continue

                # Create the vote.
                vote = Vote(author=author, post_id=post_id, type=vote_type, date=vote_date)

                yield vote

        # Insert votes in batch. Bypasses the signals!
        Vote.objects.bulk_create(vote_generator(), batch_size=1000)

        log("migrated %s votes" % Vote.objects.all().count())

    def handle(self, *args, **options):

        source = os.path.expanduser(options['dir'])

        if options['users']:
            fname = path_join(source, "users.txt")
            self.migrate_users(source, fname)

        if options['posts']:
            fname = path_join(source, "posts.txt")
            self.migrate_posts(source, fname)

        if options['votes']:
            fname = path_join(source, "votes.txt")
            self.migrate_votes(source, fname)

        if options['tags']:
            fname = path_join(source, "posts.txt")
            self.auto_tag(source, fname)

        print("*** migrated from %s" % source)
