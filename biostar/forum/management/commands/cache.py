import requests
import hjson
import time
import logging
import os
from calendar import Calendar
from itertools import chain
from urllib.parse import urljoin

from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.accounts.models import User

logger = logging.getLogger(settings.LOGGER_NAME)


def get_data(full_url):
    """
    Fetch data from remote url
    """
    while True:

        try:
            # 5 sec timeout then retry
            response = requests.get(full_url, timeout=5)
            logger.info(f"Hit remote site:{full_url}")
            break
        except Exception as exc:
            logger.error(f"{exc}...sleeping for 5 seconds then retrying.")
            time.sleep(2)

    data = {} if response.status_code == 404 else hjson.loads(response.text)
    return data, response


def cache_post(root_dir, post_list):

    api_url = "https://www.biostars.org/api/post/"

    for post_id in post_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{post_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"Post Id {post_id} does not exist.")
            continue
        # Get cache file name
        uid = data.get("id")
        fname = f"post_{uid}.hjson"

        # Write data to file
        fullpath = os.path.join(root_dir, fname)
        open(fullpath, "w").write(hjson.dumps(data))
        logger.info(f"Cached post={post_id} into {fullpath}")

    return


def cache_users(root_dir, user_list):

    api_url = "https://www.biostars.org/api/user/"
    for user_id in user_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{user_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"User Id {user_id} does not exist.")
            continue
        # Get cache file name
        uid = data.get("id")
        fname = f"user_{uid}.hjson"

        # Write data to file
        fullpath = os.path.join(root_dir, fname)
        open(fullpath, "w").write(hjson.dumps(data))
        logger.info(f"Cached user={user_id} into {fullpath}")

    return


def cache_votes(root_dir, votes_list):

    api_url = "https://www.biostars.org/api/vote/"
    for vote_id in votes_list:
        # Build url and fetch data for post
        full_url = urljoin(api_url, f"{vote_id}")
        data, response = get_data(full_url=full_url)

        # No data found for given post id
        if not data:
            logger.warning(f"Vote Id {vote_id} does not exist.")
            continue
        # Get cache file name
        uid = data.get("id")
        fname = f"vote_{uid}.hjson"

        # Write data to file
        fullpath = os.path.join(root_dir, fname)
        open(fullpath, "w").write(hjson.dumps(data))
        logger.info(f"Cached vote={vote_id} into {fullpath}")

    return


def cache_forum(root_dir, today=False, end_year=2009, start_year=2019):
    """
    Create user cache in root_dir from the start_date
    """
    # Get the remote site url
    api_url = "https://www.biostars.org/api/stats/date/"

    # Tuple with (year, month)
    months = [list(zip([year] * 12, range(1, 13))) for year in range(start_year, end_year + 1, -1)]
    # Flatten list
    months = list(chain(*months))
    # List of all days in given time period
    dates = [list(Calendar().itermonthdates(year=year, month=month)) for year, month in months]
    dates = list(chain(*dates))

    # Format the date to add extra '0' in the front.
    format_date = lambda digit: f"0{digit}" if len(f"{digit}") == 1 else f"{digit}"
    count = 0
    for date in dates:

        # Fetch the stats of the day
        month = format_date(date.month)
        day = format_date(date.day)
        daystr = f"{date.year}/{month}/{day}"
        full_url = urljoin(api_url, daystr)
        data, response = get_data(full_url=full_url)

        if not data:
            logger.error(f"No data for date={day}")
            continue

        # Create directory to store cache in
        dirname = os.path.join(root_dir, f"{date.year}", month, day)
        os.makedirs(dirname, exist_ok=True)

        # Get posts of the day
        cache_post(root_dir=dirname, post_list=data.get("new_posts"))

        # Get users of the day
        cache_users(root_dir=dirname, user_list=data.get("new_users"))

        # Get votes of the day
        cache_votes(root_dir=dirname, votes_list=data.get("new_votes"))

        # TODO: take out after testing
        count += 1
        if count == 15:
            break
    return


class Command(BaseCommand):
    help = 'Cache posts and users from remote site.'

    def add_arguments(self, parser):
        parser.add_argument("--root", help="Root directory to store into.")
        pass

    def handle(self, *args, **options):

        root = options["root"]
        root = os.path.abspath(root)
        os.makedirs(root, exist_ok=True)

        cache_forum(root_dir=root)
        pass