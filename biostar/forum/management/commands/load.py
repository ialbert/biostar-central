

import logging
import hjson
import os
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.forum.auth import create_post_from_json, create_vote
from biostar.forum.models import Post, Vote
from biostar.accounts.auth import create_user_from_json
from biostar.accounts.models import Profile


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

    vote = create_vote(author=user, post=post, vote_type=int(vote_type))

    return vote


def load_objects_from_file(root, source_file, limit, loader, is_json=False):
    """

    Parse files and pass onto the "loader"
    which is the main function used to create objects.

    """

    stream = open(source_file, 'r')
    loaded = 0

    for spec in stream:
        if loaded == limit:
            break

        if is_json:
            # Open separate json file corresponding to a single object: post/user
            json_path = os.path.join(root, spec.strip())
            json_stream = open(json_path, "r")
            info = hjson.load(json_stream)
        else:
            # Current line corresponds to a single object: vote
            # next(stream) is called to skip the header ( first line ).
            info = spec.strip() if loaded > 0 else next(stream).strip()

        obj = loader(info)
        logger.info(f"Created {obj} from {spec}")
        loaded += 1

    return


class Command(BaseCommand):
    help = 'Add existing posts/users/votes to the forum'

    def add_arguments(self, parser):
        parser.add_argument("--n",
                            type=int,
                            default=10, help="How many objects ( users,posts, votes) to load")
        parser.add_argument("--root",
                            help="Root directory with posts, votes, and users.",
                            required=True)
        parser.add_argument("--posts",
                            help="Text file with each line being a relative path to one json file ( a post )."
                            )
        parser.add_argument("--votes",
                            help="TSV with the following header: user_id post_id votes_type")
        parser.add_argument("--users",
                            help="Text file with each line being a relative path to one json file ( a user)."
                            )

    def handle(self, *args, **options):

        root = options["root"]
        nobjs = options["n"]
        posts = options["posts"]
        votes = options["votes"]
        users = options["users"]

        if users:
            users_file = os.path.join(root, users)
            load_objects_from_file(root=root, source_file=users_file, is_json=True,
                                   limit=nobjs, loader=create_user_from_json)

        if posts:
            post_file = os.path.join(root, posts)
            load_objects_from_file(root=root, source_file=post_file,is_json=True,
                                   limit=nobjs, loader=create_post_from_json)

        if votes:
            votes_file = os.path.join(root, votes)
            load_objects_from_file(root=root, source_file=votes_file,
                                   limit=nobjs, loader=vote_loader)
























