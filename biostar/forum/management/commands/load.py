

import logging

import hjson
import os
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.forum.auth import create_post_from_json, create_vote_from_json
from biostar.accounts.auth import create_user_from_json


logger = logging.getLogger(settings.LOGGER_NAME)


def load_objects_from_file(root, source_file, limit, loader):
    """

    Parse json files into dict objects and pass onto the "loader"
    which is the main function used to create objects.

    """

    stream = open(source_file, 'r')
    current = 0
    for spec in stream:
        if current == limit:
            break
        # Open the json file corresponding to a single object: post/vote/user
        json_path = os.path.join(root, spec.strip())
        json_stream = open(json_path, "r")

        info_dict = hjson.load(json_stream)
        obj = loader(info_dict)

        logger.info(f"Created {obj} id={obj.id} from {json_path}")
        json_stream.close()
        current += 1

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
                            help="Csv the following header: user_id, post_id, votes_type")
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
            load_objects_from_file(root=root, source_file=users_file,
                                   limit=nobjs, loader=create_user_from_json)

        if posts:
            post_file = os.path.join(root, posts)
            load_objects_from_file(root=root, source_file=post_file,
                                   limit=nobjs, loader=create_post_from_json)

        if votes:
            votes_file = os.path.join(root, votes)
            load_objects_from_file(root=root, source_file=votes_file,
                                   limit=nobjs, loader=create_vote_from_json)
























