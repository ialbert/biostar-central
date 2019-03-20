

import logging
import hjson
import os
import time
import requests
from urllib.parse import urljoin
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.forum.auth import create_post_from_json, preform_vote
from biostar.forum.models import Post, Vote
from biostar.accounts.auth import create_user_from_json
from biostar.accounts.models import Profile, User


logger = logging.getLogger(settings.LOGGER_NAME)


def vote_loader(single_vote):

    user_uid, post_uid, vote_type = single_vote.split("\t")

    user = Profile.objects.filter(uid=user_uid).first().user
    post = Post.objects.filter(uid=post_uid).first()
    vote = Vote.objects.filter(author=user, post=post, type=vote_type)

    if vote:
        logger.error(f"""Vote with user.uid={user_uid}, 
                    post.uid={post_uid} and vote_type={vote_type} already exists""")
        return vote.first()

    vote = preform_vote(user=user, post=post, vote_type=int(vote_type))

    return vote


def load_votes(source_file, limit=100):
    stream = open(source_file, 'r')
    loaded = 0

    for spec in stream:
        if loaded == limit:
            break
        # Current line corresponds to a single object: vote
        # next(stream) is called to skip the header ( first line ).
        info = spec.strip() if loaded > 0 else next(stream).strip()

        vote_loader(single_vote=info)

    return


def load_posts(root, source_file, limit=100):

    stream = open(source_file, 'r')
    loaded = 0

    for spec in stream:
        if loaded == limit:
            break
        json_path = os.path.join(root, spec.strip())
        json_stream = open(json_path, "r")
        info = hjson.load(json_stream)

        create_post_from_json(json_dict=info)

    return


def load_users(root, source_file, limit):

    stream = open(source_file, 'r')
    loaded = 0

    for spec in stream:
        if loaded == limit:
            break
        json_path = os.path.join(root, spec.strip())
        json_stream = open(json_path, "r")
        info = hjson.load(json_stream)

        create_user_from_json(json_dict=info)
    return


def create_post(postid):
    """
    Create post from
    """

    # Get parent post and check if it exists,
    # recursively create parent post before proceeding to current post.
    api_url = "https://www.biostars.org/api/post/"

    full_url = urljoin(api_url, f"{postid}")
    return


def posts_from_api():

    api_url = "https://www.biostars.org/api/post/"

    nposts = 400000
    #TODO: need to change listing
    for postids in range(nposts, 1, -1):

        # Get api url for user
        full_url = urljoin(api_url, f"{postids}")

        # 5 second time delay every 10 users to avoid overloading remote site.
        if postids % 10 == 0:
            time.sleep(5)

        response = requests.get(full_url)
        data = response.text

        if not data or response.status_code == 404:
            continue

        uid = data.get("id", "")
        author_uid = data.get("author_id", "")
        has_accepted = data.get("has_accepted", "")

        post = Post.objects.filter(uid=uid)
        # Check if parent exists and create that first.

        # Update an existing post
        if post:
            continue
        # Create a new post





    return




class Command(BaseCommand):
    help = 'Add existing posts/users/votes to the forum'

    def add_arguments(self, parser):
        parser.add_argument("--n",
                            type=int,
                            default=10, help="How many objects ( users,posts, votes) to load")
        parser.add_argument("--root",
                            help="Root directory with files.",
                            required=True)
        parser.add_argument("--posts",
                            help="Text file with each line being a relative path to one json file ( a post )."
                            )
        parser.add_argument("--votes",
                            help="TSV with the following header: user_id post_id votes_type")
        parser.add_argument("--users",
                            help="Text file with each line being a relative path to one json file ( a user)."
                            )
        parser.add_argument("--from_api", action="store_true",
                            help="Load posts and votes from remote site."
                            )
        parser.add_argument("--create_cache", action="store_true",
                            help="Create user,posts, and votes cache files from API"
                            )

    def handle(self, *args, **options):

        root = options["root"]
        nobjs = options["n"]
        posts = options["posts"]
        votes = options["votes"]
        users = options["users"]

        from_api = options["from_cache"]

        if from_api:
            posts_from_api()

        if users:
            users_file = os.path.join(root, users)
            load_users(root=root, source_file=users_file, limit=nobjs)

        if posts:
            post_file = os.path.join(root, posts)
            load_posts(root=root, source_file=post_file, limit=nobjs)

        if votes:
            votes_file = os.path.join(root, votes)
            load_votes(source_file=votes_file, limit=nobjs,)
























