"""
Takes as input a CVS file with two columns

Email, Name, Handler

"""
import csv
import logging
import os
import hjson
import time
import requests
from urllib.parse import urljoin

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.accounts import util
from biostar.accounts.models import User, Profile

logger = logging.getLogger("engine")


def generate_info(name):
    """
    Generate a unique email, username, and password from a name
    """
    uid = util.get_uuid(5)

    username = f"{uid}" + "_".join(name)
    email = username + "@testmail.com"
    password = f"{util.get_uuid(16)}"

    return email, username, password


def user_from_api():
    """Update or create user from remote API"""

    api_url = "https://www.biostars.org/api/user/"
    # TODO: need to change listing
    nusers = 100000

    for userids in range(nusers, 1, -1):

        # Get api url for user
        full_url = urljoin(api_url, f"{userids}")

        # 5 second time delay every 10 users to avoid overloading remote site.
        if userids % 10 == 0:
            time.sleep(5)

        response = requests.get(full_url)
        data = hjson.loads(response.text)

        # No data found for the given user id
        if not data or response.status_code == 404:
            continue

        # Get user from uid
        uid = data.get("id", "")
        user = User.objects.filter(profile__uid=uid).first()

        # Update existing user information.
        if user:
            continue

        # Create a new user with the using name and id
        name = data.get("name", "")
        email, username, password = generate_info(name=name)
        user = User.objects.create(username=username, email=email, password=password)

        # Update the profile with correct user id.
        Profile.objects.filter(user=user).update(uid=uid)

        print(userids, full_url)
        print(data)
        1 / 0

    return


class Command(BaseCommand):
    help = "Add users"

    def add_arguments(self, parser):

        parser.add_argument('fname', help="The CSV file with the users to be added. Must have headers: Name, Email")

    def handle(self, *args, **options):

        # Get the filename.
        fname = options['fname']

        if not os.path.isfile(fname):
            logger.error(f'Not a valid filename: {fname}')
            return

        stream = open(fname, 'rU')
        for row in csv.DictReader(stream):
            name = row.get("Name", '')
            email = row.get("Email", '')
            handler = row.get("Handler", '')

            # All fields must be filled in.
            if not (name and email):
                logger.error(f"Invalid row: {row}")
                continue

            # Fix capitalization.
            name = name.strip()
            email = email.strip().lower()
            handler = handler.strip()

            if User.objects.filter(email=email).exists():
                logger.info(f"Skipped creation. User with email={email} already exists.")
            else:
                username = handler or util.get_uuid(16)
                user = User.objects.create(email=email, username=username, first_name=name)
                user.set_password(settings.SECRET_KEY)
                user.save()
                logger.info(f"Created user name={name} email={email}")
