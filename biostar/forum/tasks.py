import functools
import random
from biostar.accounts.tasks import create_messages
from biostar.emailer.tasks import send_email
from django.conf import settings
import time, random
from biostar.utils.decorators import spooler, threaded
from biostar.celery import celery_task


from django.db.models import Q
#
# Do not use logging in tasks! Deadlocking may occur!
#
# https://github.com/unbit/uwsgi/issues/1369

if settings.TASKS_CELERY:
    task = celery_task
elif settings.MULTI_THREAD:
    task = threaded
else:
    task = spooler


def message(msg, level=0):
    print(f"{msg}")


@task
def classify_spam(uid):
    """
    Score the spam with a slight delay.
    """
    from biostar.forum import spam
    from biostar.forum.models import Post

    post = Post.objects.filter(uid=uid).first()

    # Give spammers the illusion of success with a slight delay
    time.sleep(1)

    # Non spam posts are left alone.
    if post.not_spam:
        return

    try:
        # Give this post a spam score and quarantine it if necessary.
        spam.score(uid=uid)

        # Add this post to the spam index.
        spam.add_spam(uid=uid)

    except Exception as exc:
        message(exc)


@task
def notify_watched_tags(uid, extra_context):
    """
    Notify users watching a given tag found in post.
    """
    from biostar.accounts.models import User
    from biostar.forum.models import Post
    from django.conf import settings

    post = Post.objects.filter(uid=uid).first()

    # Update template context with post
    extra_context.update(dict(post=post))

    users = [User.objects.filter(profile__watched__name__iexact=tag.name).distinct()
             for tag in post.root.tags.all()]

    # Flatten nested users queryset and get email.
    emails = set(u.email for o in users for u in o)

    from_email = settings.DEFAULT_NOREPLY_EMAIL

    send_email(template_name='messages/watched_tags.html',
               extra_context=extra_context,
               name=post.author.profile.name,
               recipient_list=emails,
               from_email=from_email,
               mass=True)


@task
def update_spam_index(uid):
    """
    Update spam index with this post.
    """
    from biostar.forum import spam, models

    post = models.Post.objects.filter(uid=uid).first()

    # Index posts explicitly marked as SPAM or NOT_SPAM
    # indexing SPAM increases true positives.
    # indexing NOT_SPAM decreases false positives.
    if post.spam == models.Post.DEFAULT:
        return

    # Update the spam index with most recent spam posts
    try:
        spam.add_spam(post=post)
    except Exception as exc:
        message(exc)


@task
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

# Set in the settings.
# if celery:
#     task = app.task
# elif uwsgi:
#     task = spool
# else:
#     task = threaded
#
#@spool(pass_arguments=True)
# Do this with celery.
#@shared_task
#@task
@task
def create_user_awards(user_id):

    from biostar.accounts.models import User
    from biostar.forum.models import Award, Badge, Post
    from biostar.forum.awards import ALL_AWARDS
    from biostar.forum import util

    user = User.objects.filter(id=user_id).first()
    # debugging
    # Award.objects.all().delete()

    # Collect valid targets
    valid = []

    for award in ALL_AWARDS:

        # Valid award targets the user has earned
        targets = award.validate(user)
        for target in targets:

            date = util.now()
            post = target if isinstance(target, Post) else None
            badge = Badge.objects.filter(name=award.name).first()

            # Do not award a post multiple times.
            already_awarded = Award.objects.filter(user=user, badge=badge, post=post).exists()
            if post and already_awarded:
                continue

            valid.append((user, badge, date, post))

    # Pick random awards to give to user
    random.shuffle(valid)

    valid = valid[:settings.MAX_AWARDS]

    for target in valid:
        user, badge, date, post = target

        # Create an award for each target.
        Award.objects.create(user=user, badge=badge, date=date, post=post)
        message(f"award {badge.name} created for {user.email}")


@task
def mailing_list(emails, uid, extra_context={}):
    """
    Generate notification for mailing list users.
    """
    from django.conf import settings
    from biostar.forum.models import Post

    post = Post.objects.filter(uid=uid).first()

    # Update template context with post
    extra_context.update(dict(post=post))

    # Prepare the templates and emails
    email_template = "messages/mailing_list.html"
    author = post.author.profile.name
    from_email = settings.DEFAULT_NOREPLY_EMAIL

    send_email(template_name=email_template,
               extra_context=extra_context,
               name=author,
               from_email=from_email,
               recipient_list=emails,
               mass=True)


@task
def notify_followers(sub_ids, author_id, uid, extra_context={}):
    """
    Generate notification to users subscribed to a post, excluding author, a message/email.
    """
    from biostar.forum.models import Subscription
    from biostar.accounts.models import Profile, User
    from biostar.forum.models import Post
    from django.conf import settings

    # Template used to send local messages
    local_template = "messages/subscription_message.md"

    # Template used to send emails with
    email_template = "messages/subscription_email.html"

    # Does not have subscriptions.
    if not sub_ids:
        return

    post = Post.objects.filter(uid=uid).first()
    author = User.objects.filter(id=author_id).first()
    subs = Subscription.objects.filter(uid__in=sub_ids)

    users = [sub.user for sub in subs]
    user_ids = [u.pk for u in users]

    # Update template context with post
    extra_context.update(dict(post=post))

    # Every subscribed user gets local messages with any subscription type.
    create_messages(template=local_template,
                    extra_context=extra_context,
                    user_ids=user_ids,
                    sender=author)

    # Select users with email subscriptions.
    # Exclude mailing list users to avoid duplicate emails.
    email_subs = subs.filter(type=Subscription.EMAIL_MESSAGE)
    email_subs = email_subs.exclude(user__profile__digest_prefs=Profile.ALL_MESSAGES)

    # No email subscriptions
    if not email_subs:
        return

    recipient_list = [sub.user.email for sub in email_subs]
    from_email = settings.DEFAULT_NOREPLY_EMAIL

    send_email(template_name=email_template,
               extra_context=extra_context,
               name=author.profile.name,
               from_email=from_email,
               recipient_list=recipient_list,
               mass=True)
