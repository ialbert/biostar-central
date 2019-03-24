

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
from biostar.accounts.management.commands.add_user import make_user
from biostar.forum.models import Post, Vote
from biostar.accounts.auth import create_user_from_json
from biostar.accounts.models import Profile, User


logger = logging.getLogger(settings.LOGGER_NAME)


class Bunch():
    def __init__(self, **kwargs):
        self.template = self.uid = self.text = ''
        self.status = self.tag_val = self.html = ""
        self.text = self.title = self.type = ""
        self.view_count = self.creation_date = ""
        self.lastedit_date = self.parent_id = self.root_id = ""
        self.author_id = ""
        self.__dict__.update(kwargs)


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


def get_data(full_url):
    while True:

        try:
            # 5 min timeout
            response = requests.get(full_url, timeout=300)
            data = hjson.loads(response.text)
            print("hit remote site")
            break
        except Exception as exc:
            print(f"ERROR {exc}...sleeping for 5 seconds then retrying.")
            time.sleep(5)

    return data, response


def bunch_data(data):

    bunched = Bunch(status=data.get("status_id", None), tag_val=data.get("tag_val", ""),
                    html=data.get("xhtml", ""), text=util.strip_tags(data.get("xhtml", "")),
                    title=data.get("title", ""), type=data.get("type_id", None), view_count=data.get("view_count", 0),
                    creation_date=data.get("creation_date"), lastedit_date=data.get("lastedit_date"),
                    parent_id=data.get("parent_id"), root_id=data.get("root_id"), uid=data.get("id"),
                    author_id=data.get("author_id", ""))

    return bunched


def create_parent(parent_id):

    print(f"CREATING PARENT SEPERATLY: {parent_id}")
    api_url = "https://www.biostars.org/api/post/"
    full_url = urljoin(api_url, f"{parent_id}")
    print(parent_id)
    data, response = get_data(full_url=full_url)

    if not data or response.status_code == 404:
        print(f"postid {parent_id} does not exist.")
        return

    data = bunch_data(data=data)
    parent = Post.objects.filter(uid=parent_id).first()
    post = Post.objects.filter(uid=data.uid).first()
    user = User.objects.filter(profile__uid=data.author_id).first()

    # Skip posts without authors.
    if not user:
        print(f"user with {data.author_id} does not exist.")
        make_user(userid=data.author_id)
        return
    if not post:
        post = Post.objects.create(tag_val=data.tag_val, uid=data.uid, title=data.title, content=data.text,
                                   type=data.type, html=data.html, view_count=data.view_count,
                                   creation_date=data.creation_date, author=user, lastedit_date=data.lastedit_date,
                                   status=data.status, parent=parent)
    return post


def make_post(postid):
    """
    Create post from
    """

    api_url = "https://www.biostars.org/api/post/"
    full_url = urljoin(api_url, f"{postid}")
    print(postid)
    data, response = get_data(full_url=full_url)

    if not data or response.status_code == 404:
        print(f"postid {postid} does not exist.")
        return

    data = bunch_data(data=data)
    parent = Post.objects.filter(uid=data.parent_id).first()
    root = Post.objects.filter(uid=data.root_id).first()
    post = Post.objects.filter(uid=data.uid).first()

    user = User.objects.filter(profile__uid=data.author_id).first()

    # Skip posts without authors.
    if not user:
        print(f"user with {data.author_id} does not exist.")
        make_user(userid=data.author_id)
        return

    # Recursively create parent post before proceeding to current post.
    if not root and (data.root_id != data.uid):
        print(f"ROOT for {postid} does not exist. CREATING THE ROOT {data.root_id}")
        time.sleep(.5)
        make_post(postid=data.root_id)

    if not parent and data.parent_id != data.uid:
        print(f"PARENT for {postid} does not exist. CREATING THE PARENT {data.parent_id}")
        time.sleep(.5)
        make_post(postid=data.parent_id)

    # Update an existing post
    if post:
        print(f"Updating existing post {postid}")
        Post.objects.filter(uid=data.uid).update(tag_val=data.tag_val, title=data.title, content=data.text, type=data.type,
                                                 html=data.html, view_count=data.view_count, creation_date=data.creation_date,
                                                 author=user, lastedit_date=data.lastedit_date, status=data.status,
                                                 parent=parent, root=root)
    else:
        # Create a new post
        print(f"Creating new post {postid}")
        if data.type in (Post.COMMENT, Post.ANSWER) and not (parent or root):
            create_parent(parent_id=data.parent_id)
        else:
            Post.objects.create(tag_val=data.tag_val, uid=data.uid, title=data.title, content=data.text, type=data.type,
                                html=data.html, view_count=data.view_count, creation_date=data.creation_date, author=user,
                                lastedit_date=data.lastedit_date, status=data.status, parent=parent)

    return


def posts_from_api():

    #nposts = 364053
    nposts = 339359

    #TODO: need to change listing
    for postid in range(nposts, 0, -1):

        # 5 second time delay every 30 posts to avoid overloading remote site.
        if postid % 30 == 0:
            print("Entering 5s timedelay")
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
                            help="Root directory with files.")
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

        from_api = options["from_api"]

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
























