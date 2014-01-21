from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from django.utils.timezone import utc
from django.conf import settings
from django.contrib.auth import get_user_model
import os, csv, datetime
from django.utils.dateparse import parse_datetime
from django.utils import timezone


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

MIGRATE_DIR = os.environ["BIOSTAR_MIGRATE_DIR"]
MIGRATE_DIR = os.path.expanduser(MIGRATE_DIR)


def get(data, attr):
    return data.get(attr, '').strip()


def localize_time(text):
    naive = parse_datetime(text)
    local = timezone.make_aware(naive, timezone=utc)
    return local


class Command(BaseCommand):
    help = 'migrate data from Biostar 1.*'

    option_list = BaseCommand.option_list + (
        make_option("-u", '--user', action="store_true", dest='user', default=False, help='user data import file'),
    )

    def migrate_users(self, fname):

        #User.objects.all().delete()

        self.stdout.write("migrating user data from %s" % fname)
        stream = csv.DictReader(file(fname), delimiter=b'\t')

        seen = set()

        for row in stream:
            uid = get(row, 'id')

            # The first user is the default admin.
            if int(uid) == 1:
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
            prof.about_me = file(path_join(MIGRATE_DIR, 'about_me', uid)).read()

            prof.save()

            self.stdout.write("migrated %s" % user)

            #print(row)

        for user in User.objects.all():
            print (user)


    def handle(self, *args, **options):
        user_file = options['user']
        if user_file:
            user_file = path_join(MIGRATE_DIR, "users.txt")
            self.migrate_users(user_file)

