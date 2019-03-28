
import logging
import  sqlite3
import os
from django.core.management.base import BaseCommand
from biostar.accounts.models import Profile
from django.contrib.auth import get_user_model

User = get_user_model()
logger = logging.getLogger("engine")


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.uid = self.text = ''
        self.user = self.stream = None
        self.__dict__.update(kwargs)


def map_row(row, colnames):
    """
    Pair items in row with corresponding column name.
    Return dictionary keyed by column names.
    """
    mapped = {col_name: row for row, col_name in zip(row, colnames)}

    return mapped


def parse_user(row, dbfile):
    """
    Parse user info from row, return a dict keyed by colum name.
    """
    # Connect the database and get relevant user tables
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    users_table = "users_user"
    profile_table = "users_profile"

    # Get user table to extract column names from
    cursor.execute(f"SELECT * FROM {users_table}")
    user_colnames = [col[0] for col in cursor.description]

    # Get user info in as dictionary
    user = map_row(row=row, colnames=user_colnames)

    # Get profile for specific user from profile table
    cursor.execute(f"SELECT * FROM {profile_table} WHERE user_id={user['id']}")
    # Pick the first profile returned for the user.id
    profile = cursor.fetchall()[0]

    # Get column names for profile table
    profile_colnames = [col[0] for col in cursor.description]
    # Get profile info as dictionary
    user_profile = map_row(row=profile, colnames=profile_colnames)

    # Update user dict with profile info.
    user.update(user_profile)

    return user


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        parser.add_argument('--database', help="Source database to copy from.", required=True)

        pass

    def handle(self, *args, **options):

        # Get users from default database
        # TODO: use sqlalchemy for postgres
        dbfile = options["database"]
        dbfile = os.path.abspath(dbfile)

        # Connect the database
        conn = sqlite3.connect(dbfile)
        cursor = conn.cursor()

        # Get the relevant user table
        users_table = "users_user"

        # Get the user rows to iterate over.
        cursor.execute(f"SELECT * FROM {users_table}")

        for row in cursor:
            # Get the user info from row
            user = parse_user(row=row, dbfile=dbfile)

            new_user = User.objects.filter(email=user["email"]).first()
            # Skip when users that already exists.
            if new_user:
                continue

            # Create user in target database.
            username = user.get("username") or f"{user['name']}{user['id']}"
            new_user = User.objects.create(username=username, email=user["email"],
                                           password=user["password"], is_active=user["is_active"],
                                           is_superuser=user["is_admin"], is_staff=user["is_staff"])

            # Update user profile
            Profile.objects.filter(user=new_user).update(uid=user["id"], name=user["name"],
                              role=user["type"], last_login=user["last_login"],
                              date_joined=user["date_joined"], location=user["location"],
                              website=user["website"], scholar=user["scholar"], text=user["info"],
                              score=user["score"], twitter=user["twitter_id"], my_tags=user["my_tags"],
                              digest_prefs=user["digest_prefs"], new_messages=user["new_messages"],)
        return
