from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from django.utils.timezone import utc
from django.conf import settings
from django.contrib.auth import get_user_model
import os, csv, datetime
from django.utils.dateparse import parse_datetime
from django.utils import timezone, encoding
from biostar.apps.posts.models import Post, PostBody

def path_join(*args):
    return os.path.abspath(os.path.join(*args))

# Obtaine the user model
User = get_user_model()

USER_TYPE_MAP = {
    "New": User.NEW, "Member": User.MEMBER,
    "Moderator": User.MODERATOR, "Administrator": User.ADMIN
}

USER_STATUS_MAP = {
    "Active": User.ACTIVE, "Suspended": User.SUSPENDED, "Banned": User.BANNED,
}

POST_TYPE_MAP = {
    "Question": Post.QUESTION, "Answer": Post.ANSWER,
    "Comment": Post.COMMENT, "Job": Post.JOB,
}

MIGRATE_DIR = os.environ["BIOSTAR_MIGRATE_DIR"]
MIGRATE_DIR = os.path.expanduser(MIGRATE_DIR)


def get(data, attr):
    value = data.get(attr, '').strip()
    value = encoding.smart_unicode(value)
    return value

def localize_time(text):
    naive = parse_datetime(text)
    local = timezone.make_aware(naive, timezone=utc)
    return local


class Command(BaseCommand):
    help = 'migrate data from Biostar 1.*'

    option_list = BaseCommand.option_list + (
        make_option("-u", '--user', action="store_true", dest='user', default=False, help='import users'),
        make_option("-p", '--post', action="store_true", dest='post', default=False, help='import posts'),
    )

    def migrate_posts(self, fname):

        Post.objects.all().delete()

        user_map = dict((u.id, u) for u in User.objects.all())

        self.stdout.write("migrating posts data from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        for row in stream:
            uid = get(row, 'id')
            title = get(row, 'title')

            author_id = get(row, 'author_id')
            author_id = int(author_id)
            author = user_map.get(author_id)

            if not author:
                print ("*** author %s not found for post %s" % (author_id, uid))
                continue

            post_type = get(row, 'post_type')
            post_type = POST_TYPE_MAP.get(post_type, Post.FORUM)
            post_status = Post.OPEN if get(row, 'post_status') == "Open" else Post.CLOSED


            post = Post(id=uid, title=title, author=author, lastedit_user=author)
            post.status = post_status
            post.type = post_type
            post.creation_date = localize_time(get(row, 'creation_date'))
            post.lastedit_date = localize_time(get(row, 'lastedit_date'))

            post.save()

            self.stdout.write("migrated %s" % post)

    def migrate_users(self, fname):

        #User.objects.all().delete()

        self.stdout.write("migrating user data from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        seen = set()

        for row in stream:
            uid = int(get(row, 'id'))

            # The first user is the default admin.
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

            self.stdout.write("migrated %s" % user)

            #print(row)

        for user in User.objects.all():
            print (user)


    def handle(self, *args, **options):
        user_file = options['user']
        post_file = options['post']

        if user_file:
            fname = path_join(MIGRATE_DIR, "users.txt")
            self.migrate_users(fname)

        if post_file:
            fname = path_join(MIGRATE_DIR, "posts.txt")
            self.migrate_posts(fname)
