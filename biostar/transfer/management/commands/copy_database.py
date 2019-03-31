
import logging
import  sqlite3
import os
from django.core.management.base import BaseCommand
from biostar.accounts.models import User
from biostar.transfer.models import UsersUser, UsersProfile
from biostar.forum import util
logger = logging.getLogger("engine")


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.uid = self.text = ''
        self.user = self.stream = None
        self.__dict__.update(kwargs)


def copy_users():

    to_copy = UsersUser.objects.using("biostar2").all()

    for user in to_copy:
        # See if user already exists.
        new_user = User.objects.filter(email=user.email).first()

        print(user, new_user, user.profile)

        if new_user:
            continue

        new_user = User.objects.create(email=user.email, password=user.password).first()

        # Update user profile



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
