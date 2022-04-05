#
# Starter script template for writing programs that interact with the django site
#

import logging, plac
from itertools import *
from biostar.tools.config import *


logger = logging.getLogger("engine")

@plac.opt("limit", "limit ", abbrev="L")
@plac.opt("days", "how far to go back in time", abbrev="d")
@plac.flg("apply", "applies the action", abbrev='A')
@plac.flg("verbose", "show debug messages")
def main(apply=False, days=7, limit=10, verbose=False):

    # Increase verbosity
    if verbose:
        logger.setLevel(logging.DEBUG)

    # Get the first admin.
    admin = get_first_admin()

    # Get the recently created spam
    recent = time_ago(days=days)
    posts = Post.objects.filter(spam=Post.SPAM, creation_date__gt=recent).select_related("author__profile").order_by("-pk")

    # Apply the limit.
    posts = islice(posts, limit)

    # Find unique users that posted spam.
    users = {}
    for post in posts:
        logger.info(f"target: user={post.author.profile.name}, title={post.title}")
        users[post.author.id] = post.author

    # Various counter.
    user_count = delete_count = 0

    # Additional sanity checks on users
    for user in users.values():

        # Post query.
        posts = Post.objects.filter(author=user)

        # How posts for the user.
        ham_count = posts.exclude(spam=Post.SPAM).count()

        # How many spam posts.
        spam_count = posts.filter(spam=Post.SPAM).count()

        # More spam than ham
        if spam_count > ham_count:
            user_count += 1
            delete_count += spam_count
            text = f"spammer: user={user.profile.name}, spam_count={spam_count}"
            logger.info(text)
            if apply:
                posts.filter(spam=Post.SPAM).delete()

    if apply and delete_count:
        msg = f"spam cleanup, removed {user_count} spammers and {delete_count} spam posts"
        auth.db_logger(user=admin, text=msg)


if __name__ == "__main__":
    plac.call(main)
