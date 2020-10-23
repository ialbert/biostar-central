
from biostar.accounts.tasks import create_messages
from biostar.emailer.tasks import send_email
import time
from biostar.utils.decorators import spool, timer


from django.db.models import Q
#
# Do not use logging in tasks! Deadlocking may occur!
#
# https://github.com/unbit/uwsgi/issues/1369


def message(msg, level=0):
    print(f"{msg}")


@spool(pass_arguments=True)
def spam_scoring(post):
    """
    Score the spam with a slight delay.
    """
    from biostar.forum import spam

    # Give spammers the illusion of success with a slight delay
    time.sleep(1)

    try:
        # Give this post a spam score and quarantine it if necessary.
        spam.score(post=post)
    except Exception as exc:
        message(exc)


@spool(pass_arguments=True)
def notify_watched_tags(post):
    """
    Create a subscription to a post when
    """
    from biostar.accounts.models import User
    from biostar.forum import auth

    # Skip non top level posts.
    if not post.is_toplevel:
        return

    # Iterate over tags and get users that are watching them
    # Exclude mailing-list mode users.
    users = [User.objects.filter(profile__watched_tags__contains=tag) for tag in post.tags.all()]

    # Flatten nested iterable.
    users = set(u for qs in users for u in qs)

    # Subscribe users to this post and send them emails.
    for user in users:
        auth.create_subscription(post=post.root, user=user)

    return


@spool(pass_arguments=True)
def update_spam_index(post):
    """
    Update spam index with this post.
    """
    from biostar.forum import spam

    # Index posts explicitly marked as SPAM or NOT_SPAM
    # indexing SPAM increases true positives.
    # indexing NOT_SPAM decreases false positives.
    if not (post.is_spam or post.not_spam):
        return

    # Update the spam index with most recent spam posts
    try:
        spam.add_spam(post=post)
    except Exception as exc:
        message(exc)


@spool(pass_arguments=True)
def created_post(pid):
    message(f"Created post={pid}")
    pass


#
# This timer leads to problems as described in
#
# https://github.com/unbit/uwsgi/issues/1369
#

# #@timer(secs=180)
# def update_index(*args):
#     """
#     Index 1000 posts every 3 minutes
#     """
#     from biostar.forum.models import Post
#     from biostar.forum import search
#     from django.conf import settings
#
#     # Get un-indexed posts
#     posts = Post.objects.filter(indexed=False)[:settings.BATCH_INDEXING_SIZE]
#
#     # Nothing to be done.
#     if not posts:
#         message("No new posts found")
#         return
#
#     message(f"Indexing {len(posts)} posts.")
#
#     # Update indexed field on posts.
#     Post.objects.filter(id__in=posts.values('id')).update(indexed=True)
#
#     try:
#         search.index_posts(posts=posts)
#         message(f"Updated search index with {len(posts)} posts.")
#     except Exception as exc:
#         message(f'Error updating index: {exc}')
#         Post.objects.filter(id__in=posts.values('id')).update(indexed=False)
#
#     return

@spool(pass_arguments=True)
def create_user_awards(user_id):

    from biostar.accounts.models import User
    from biostar.forum.models import Award, Badge, Post
    from biostar.forum.awards import ALL_AWARDS

    user = User.objects.filter(id=user_id).first()

    # debugging
    # Award.objects.all().delete()

    for award in ALL_AWARDS:
        # Valid award targets the user has earned
        targets = award.validate(user)

        for target in targets:
            date = user.profile.last_login
            post = target if isinstance(target, Post) else None
            badge = Badge.objects.filter(name=award.name).first()

            # Do not award a post multiple times.
            already_awarded = Award.objects.filter(user=user, badge=badge, post=post).exists()
            if post and already_awarded:
                continue

            # Create an award for each target.
            Award.objects.create(user=user, badge=badge, date=date, post=post)

            message("award %s created for %s" % (badge.name, user.email))


@spool(pass_arguments=True)
def mailing_list(users, extra_context={}):
    """
    Generate notification for mailing list users.
    """

    # Prepare the templates and emails
    email_template = "messages/mailing_list.html"
    emails = [user.email for user in users]
    send_email(template_name=email_template, extra_context=extra_context, recipient_list=emails)


@spool(pass_arguments=True)
def notify_followers(subs, author, extra_context={}):
    """
    Generate notification to users subscribed to a post, excluding author, a message/email.
    """
    from biostar.forum.models import Subscription

    # Template used to send local messages
    local_template = "messages/subscription_message.md"

    # Template used to send emails with
    email_template = "messages/subscription_email.html"

    # Does not have subscriptions.
    if not subs:
        return

    # Select users that should be notified.
    # Exclude mailing list users to avoid duplicate emails.
    users = [sub.user for sub in subs if not sub.user.profile.mailing_list]

    # Every subscribed user gets local messages with any subscription type.
    create_messages(template=local_template, extra_context=extra_context, rec_list=users, sender=author)

    # Select users with email subscriptions.
    email_subs = subs.filter(type=Subscription.EMAIL_MESSAGE)

    # No email subscriptions
    if not email_subs:
        return

    recipient_list = [sub.user.email for sub in email_subs]

    send_email(template_name=email_template, extra_context=extra_context, recipient_list=recipient_list)
