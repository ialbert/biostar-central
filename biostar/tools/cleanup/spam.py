#
# Starter script template for writing programs that interact with the django site
#

import logging, plac
from itertools import *
from biostar.tools.config import *

logger = logging.getLogger("biostar")

@plac.opt("limit", "limit ", abbrev="L")
@plac.opt("days", "deletes the users", abbrev="d")
@plac.flg("apply", "applies the action", abbrev='D')
@plac.flg("verbose", "show debug messages")
def main(apply=False, days=2, limit=10, verbose=False):
    # Increase verbosity
    if verbose:
        logger.setLevel(logging.DEBUG)

    admin = get_first_admin()

    # Get the spam posts
    posts = Post.objects.filter(spam=Post.SPAM, author__profile__state=Profile.NEW).order_by("-pk")

    # Apply the limit.
    posts = islice(posts, limit)

    # Find unique users that posted spam.
    users = {}
    for post in posts:
        logger.debug(f"spam: {post.title}")
        users[post.author.id] = post.author

    # Various counter.
    u_count = s_count = 0

    # Additional sanity checks on users
    for user in users.values():

        # How many posts did the user make.
        posts = get_posts(user)
        post_count = posts.exclude(spam=Post.SPAM).count()

        # How many spam posts did the user make.
        spam_count = posts.filter(spam=Post.SPAM).count()

        if spam_count > post_count:
            u_count += 1
            s_count += spam_count
            text = f"user={user.profile.name} spam_count={spam_count}"
            logger.info(text)
            if apply:
                user.delete()

    if apply and s_count:
        msg = f"spam cleanup, removed {u_count} spammers and {s_count} spam posts"
        auth.db_logger(user=admin, text=msg)


if __name__ == "__main__":
    plac.call(main)
