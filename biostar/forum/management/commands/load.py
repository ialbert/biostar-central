

import logging
import hjson
import os
import time
import requests
from urllib.parse import urljoin
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.forum.auth import create_post_from_json, preform_vote
from biostar.forum import util
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


def make_post(postid):
    """
    Create post from
    """

    # Get parent post and check if it exists,
    # recursively create parent post before proceeding to current post.
    api_url = "https://www.biostars.org/api/post/"
    full_url = urljoin(api_url, f"{postid}")
    response = requests.get(full_url)
    data = hjson.loads(response.text)

    if not data or response.status_code == 404:
        return

    status = data.get("status_id", None)
    tag_val = data.get("tag_val", "")
    html = data.get("xhtml", "")
    text = util.strip_tags(html)
    title = data.get("title", "")
    post_type = data.get("type_id", None)
    view_count = data.get("view_count", 0)
    creation_date = data.get("creation_date")
    lastedit_date = data.get("lastedit_date")

    parent_id = data.get("parent_id")
    root_id = data.get("root_id")
    uid = data.get("id")
    parent = Post.objects.filter(uid=parent_id)
    root = Post.objects.filter(uid=root_id)
    post = Post.objects.filter(uid=uid)

    author_id = data.get("author_id", "")
    user = User.objects.filter(profile__uid=author_id).first()

    # Skip posts without authors.
    if not user:
        return

    # Recursively create parent/root if it does not already exist.
    if not root and root_id != uid:
        time.sleep(5)
        make_post(postid=root_id)

    if not parent and parent_id != uid:
        time.sleep(5)
        make_post(postid=parent_id)

    # Update an existing post
    if post:
        Post.objects.filter(uid=uid).update(tag_val=tag_val, title=title, text=text, type=post_type,
                                            view_count=view_count, creation_date=creation_date,
                                            lastedit_date=lastedit_date, status=status, parent=parent,
                                            root=root)
    else:
        # Create a new post
        Post.objects.create(tag_val=tag_val, uid=uid, title=title, text=text, type=post_type,
                            view_count=view_count, creation_date=creation_date,
                            lastedit_date=lastedit_date, status=status, parent=parent,
                            root=root)

    return


def posts_from_api():

    nposts = 400000
    #TODO: need to change listing
    for postid in range(nposts, 1, -1):

        # 5 second time delay every 10 users to avoid overloading remote site.
        if postid % 10 == 0:
            time.sleep(5)

        make_post(postid=postid)

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

    def handle(self, *args, **options):

        root = options["root"]
        nobjs = options["n"]
        posts = options["posts"]
        votes = options["votes"]
        users = options["users"]

        from_api = options["from_cache"]

        if from_api:
            posts_from_api()
            return

        if users:
            users_file = os.path.join(root, users)
            load_users(root=root, source_file=users_file, limit=nobjs)

        if posts:
            post_file = os.path.join(root, posts)
            load_posts(root=root, source_file=post_file, limit=nobjs)

        if votes:
            votes_file = os.path.join(root, votes)
            load_votes(source_file=votes_file, limit=nobjs,)
























