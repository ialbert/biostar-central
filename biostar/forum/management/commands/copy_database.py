
import logging
import  sqlite3
import os
from django.core.management.base import BaseCommand
from biostar.accounts.models import Profile, User
from biostar.forum.models import Post
from biostar.forum import util
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
    users = conn.cursor()

    # Get user table to extract column names from
    users_table = "users_user"
    profile_table = "users_profile"
    users.execute(f"SELECT * FROM {users_table}")
    user_colnames = [col[0] for col in users.description]

    # Get user info in as dictionary
    user = map_row(row=row, colnames=user_colnames)

    # Get profile for specific user from profile table
    users.execute(f"SELECT * FROM {profile_table} WHERE user_id={user['id']}")
    # Pick the first profile returned for the user.id
    profile = users.fetchone()
    # Get column names for profile table
    profile_colnames = [col[0] for col in users.description]
    # Get profile info as dictionary
    user_profile = map_row(row=profile, colnames=profile_colnames)
    # Update user dict with profile info.
    user.update(user_profile)
    conn.close()
    return user


def copy_users(dbfile):

    conn = sqlite3.connect(dbfile)
    users = conn.cursor()
    # Get the user table
    users_table = "users_user"
    # Get the users to iterate over.
    users.execute(f"SELECT * FROM {users_table}")

    for row in users:
        # Parse user information into dictionary
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
        Profile.objects.filter(user=new_user).update(uid=user["id"], name=user["name"], role=user["type"],
                                                     last_login=user["last_login"], date_joined=user["date_joined"],
                                                     location=user["location"], website=user["website"],
                                                     scholar=user["scholar"],  score=user["score"], text=user["info"],
                                                     twitter=user["twitter_id"], digest_prefs=user["digest_prefs"],
                                                     new_messages=user["new_messages"], my_tags=user["my_tags"],)
        print(new_user.profile.name, new_user.profile.uid)
    conn.close()
    return


def copy_posts(dbfile):

    conn = sqlite3.connect(dbfile)
    posts = conn.cursor()
    # Get the posts table
    post_table = "posts_post"
    # Get the users to iterate over.
    posts.execute(f"SELECT * FROM {post_table}")

    for row in posts:
        colnames = [col[0] for col in posts.description]

        # Get post info in as dictionary
        post = map_row(row=row, colnames=colnames)
        new_post = Post.objects.get_all(uid=post["id"]).first()

        # Skip if posts exists.
        if new_post:
            continue

        # Get the foreign key relationship
        content = util.strip_tags(post["content"])
        new_post = Post.objects.create(uid=post["id"], html=post["html"], type=post["type"],
                                       subs_count=post["subs_count"], lastedit_user_id=post["lastedit_user_id"],
                                       author_id=post["author_id"], root_id=post["root_id"], status=post["status"],
                                       parent_id=post["parent_id"], lastedit_date=post["lastedit_date"],
                                       content=content, comment_count=post["comment_count"],
                                       has_accepted=post["has_accepted"], title=post["title"],
                                       thread_score=post["thread_score"], vote_count=post["vote_count"],
                                       creation_date=post["creation_date"], tag_val=post["tag_val"],
                                       reply_count=post["reply_count"], book_count=post["book_count"],
                                       view_count=post["view_count"])
        print(new_post.title, "CREATED")
    conn.close()
    return


def copy_votes(dbfile):

    return


def copy_subs(dbfile):
    return


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

        # Copy users, posts, votes, then subscriptions in order.
        #copy_users(dbfile=dbfile)
        copy_posts(dbfile=dbfile)
        copy_votes(dbfile=dbfile)
        copy_subs(dbfile=dbfile)

        return
