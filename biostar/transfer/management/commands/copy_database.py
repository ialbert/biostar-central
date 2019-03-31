
import logging
import  sqlite3
import os
from django.core.management.base import BaseCommand
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost
from biostar.forum import util
logger = logging.getLogger("engine")


def copy_users():

    source = UsersUser.objects.all()

    for user in source:
        # See if user already exists.
        new_user = User.objects.filter(email=user.email).first()
        if new_user:
            continue

        username = f"{user.name}{user.id}" if User.objects.filter(username=user.name).exists() else user.name
        # Create user
        new_user = User.objects.create(username=username, email=user.email,
                                       password=user.password, is_active=user.is_active,
                                       is_superuser=user.is_admin, is_staff=user.is_staff)

        # Update profile
        Profile.objects.filter(user=new_user).update(uid=user.id, name=user.name,
                          role=user.type, last_login=user.last_login,
                          date_joined=user.profile.date_joined, location=user.profile.location,
                          website=user.profile.website, scholar=user.profile.scholar, text=user.profile.info,
                          score=user.score, twitter=user.profile.twitter_id, my_tags=user.profile.my_tags,
                          digest_prefs=user.profile.digest_prefs, new_messages=user.new_messages)
        logger.info(f"Created user email={new_user.email}")


def copy_posts():
    """
    Bulk create posts
    """

    source = PostsPost.objects.all()

    def generate():
        return


    print(generate())
    return


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        # Get users from default database

        # Copy users, posts, votes, then subscriptions in order.
        copy_users()
        #copy_posts()
        #copy_votes()


        return
