

import logging

import hjson
import os
from django.conf import settings
from django.core.management.base import BaseCommand


from biostar.accounts import auth as accounts_auth


logger = logging.getLogger(settings.LOGGER_NAME)



def load_posts(posts_file, limit):
    "Load posts made by users"




    return





def load_votes(votes_file, limit):
    ""


    return




def load_users(root, users_file, limit, pwd=None):
    "Recreate users from a users_file and gives random uid as password"


    stream = open(users_file, 'r')


    for user_spec in stream:

        # open the json file corresponding to a user
        json_path = os.path.join(root, user_spec)
        json_stream = open(json_path, "r")

        user_info = hjson.load(json_stream)

        user = accounts_auth.create_user_from_json(user_info)

        logger.info(f"Added name={user.profile.name} id={user.id} uid={user.profile.uid} from {json_path}")
        json_stream.close()

        1/0

    return


class Command(BaseCommand):
    help = 'Add existing posts/users/votes to the forum'

    def add_arguments(self, parser):
        parser.add_argument("--n",
                            type=int,
                            default=10, help="How many objects ( users,posts, votes) to load")
        parser.add_argument("--root",
                            help="Root directory with with posts, votes, and users.",
                            required=True)
        parser.add_argument("--posts",
                            help="Text file with each line being the root_dir/post-id." )
        parser.add_argument("--votes",
                            help="Text file following header: user_id, post_id, votes_type")
        parser.add_argument("--users",
                            help="Text file with each line being the root_dir/post-id.")

    def handle(self, *args, **options):

        # Collect the parameters.

        root = options["root"]
        nobjs = options["n"]
        posts = options["posts"]
        votes = options["votes"]
        users = options["users"]

        if users:
            users_file = os.path.join(root, users)
            load_users(users_file=users_file, limit=nobjs, root=root)

        if posts:
            post_file = os.path.join(root, posts)
            load_posts(posts_file=post_file, limit=nobjs)

        if votes:
            votes_file = os.path.join(root, votes)
            load_votes(votes_file=votes_file, limit=nobjs)






















