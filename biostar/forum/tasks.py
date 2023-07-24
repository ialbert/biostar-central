import functools
import random, logging, os
from biostar.accounts.tasks import create_messages
from biostar.emailer.tasks import send_email
from django.conf import settings
import time, random
from biostar.utils.decorators import task, timer

from django.db.models import Q

logger = logging.getLogger("engine")


#
# Do not use logging in tasks! Deadlocking may occur!
#
# https://github.com/unbit/uwsgi/issues/1369

def message(msg, level=0):
    print(f"{msg}")


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
    if emails:
        send_email(template_name='messages/watched_tags.html',
                   extra_context=extra_context,
                   name=post.author.profile.name,
                   recipient_list=emails,
                   from_email=from_email,
                   mass=True)


@task
def created_post(pid):
    message(f"Created post={pid}")
    pass


# @timer(2)
# def inner_timer(*args, **kwargs):
#     print("TIMERRRRR " *10)
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
# @spool(pass_arguments=True)
# Do this with celery.
# @shared_task
# @task
@task
def create_user_awards(user_id, limit=None):
    from biostar.accounts.models import User
    from biostar.forum.models import Award, Badge, Post
    from biostar.forum.awards import ALL_AWARDS
    from biostar.forum import util, auth
    from django.conf import settings

    limit = limit or settings.MAX_AWARDS

    user = User.objects.filter(id=user_id).first()
    # debugging
    # Award.objects.all().delete()

    # Collect valid targets
    valid = auth.valid_awards(user=user)

    # Pick random awards to give to user
    random.shuffle(valid)

    valid = valid[:limit]

    for target in valid:
        user, badge, date, post = target

        # Set the award date to the post edit date
        date = post.lastedit_date if post else date

        # Create an award for each target.
        Award.objects.create(user=user, badge=badge, date=date, post=post)

        message(f"award {badge.name} created for {user.email}")


def batch_create_awards(limit=100):
    from biostar.accounts.models import User
    from biostar.forum import auth, models, util

    # Randomly order awards
    users = User.objects.order_by('?')[:limit]

    # Aggregate target awards into flat list
    targets = []
    for u in users:
        valid = auth.valid_awards(user=u)
        targets.extend(valid)

    def batch():
        for target in targets:
            if not target:
                continue
            user, badge, date, post = target
            award = models.Award(user=user, badge=badge, date=date, post=post)
            logger.debug(f"awarded {award.badge.name} to {user}")

            yield award

    models.Award.objects.bulk_create(objs=batch(), batch_size=limit)
    logger.info(f"{len(targets)} awards given to {len(users)} users")


def high_trust(user, minscore=50):
    """
    Conditions for trusting a user
    """
    prof = user.profile
    cond = prof.is_moderator or prof.score >= minscore
    return cond


def low_trust(user, minscore=50):
    return not high_trust(user, minscore=minscore)


@task
def set_link_title(pk):
    """
    Sets the title of a shared link.
    """
    import requests
    from bs4 import BeautifulSoup
    from biostar.forum.models import SharedLink
    link = SharedLink.objects.filter(pk=pk).first()

    logger.info(f"getting link title for {link.url}")

    try:
        # Fetch the page.
        resp = requests.get(link.url)

        # Parse the content
        soup = BeautifulSoup(resp.text, 'html.parser')

        # Set the title
        for elem in soup.find_all('title'):
            title = elem.get_text()
            if title:
                title = title.strip()
                SharedLink.objects.filter(pk=pk).update(title=title)
                break

    except Exception as exc:
        logger.warning(exc)



@task
def spam_check(uid):
    from biostar.forum.models import Post, Log, delete_post_cache
    from biostar.accounts.models import User, Profile
    from biostar.forum.auth import db_logger

    post = Post.objects.filter(uid=uid).first()
    author = post.author

    if not settings.CLASSIFY_SPAM:
        return

    # Automated spam disabled in for trusted user
    if author.profile.trusted or author.profile.score > 50:
        return

    # Classify spam only if we have not done it yet.
    if post.spam != Post.DEFAULT:
        return

    # Drop the cache for the post.
    delete_post_cache(post)

    try:
        from biostar.utils import spamlib

        if not os.path.isfile(settings.SPAM_MODEL):
            spamlib.build_model(fname=settings.SPAM_DATA, model=settings.SPAM_MODEL)

        # Short posts do not get classified too many false positives
        if len(post.content) < 150:
            return

        # Classify the content.
        flag = spamlib.classify_content(post.content, model=settings.SPAM_MODEL)

        # Another process may have already classified it as spam.
        check = Post.objects.filter(uid=post.uid).first()
        if check and check.spam == Post.SPAM:
            return

        ## Links in title usually mean spam.
        spam_words = ["http://", "https://" ]
        for word in spam_words:
            flag = flag or (word in post.title.lower())
        
        spam_words2 = ["cialis", "viagra" ]
        for word in spam_words2:
            flag = flag or (word in post.title.lower() + post.content.lower())

        # Handle the spam.
        if flag:

            Post.objects.filter(uid=post.uid).update(spam=Post.SPAM, status=Post.CLOSED)

            # Get the first admin.
            user = User.objects.filter(is_superuser=True).order_by("pk").first()

            create_messages(template="messages/spam-detected.md",
                            extra_context=dict(post=post),
                            user_ids=[post.author.id])

            spam_count = Post.objects.filter(spam=Post.SPAM, author=author).count()

            db_logger(user=user, action=Log.CLASSIFY, target=post.author, text=f"classified the post as spam",
                      post=post)

            if spam_count > 1 and low_trust(post.author):
                # Suspend the user
                Profile.objects.filter(user=author).update(state=Profile.SUSPENDED)
                db_logger(user=user, action=Log.MODERATE, text=f"suspended", target=post.author)

    except Exception as exc:
        print(exc)
        logger.error(exc)

    return False


@task
def herald_emails(uid):
    """
    Send emails to herald subscribers
    """
    from biostar.emailer.models import EmailSubscription, EmailGroup
    from biostar.forum.models import Post
    post = Post.objects.filter(uid=uid).first()
    group = EmailGroup.objects.filter(uid='herald').first()
    # Get active subscriptions to herald.
    subs = EmailSubscription.objects.filter(group=group, state=EmailSubscription.ACTIVE)

    if not subs:
        return

    emails = subs.values_list('email', flat=True)
    context = dict(post=post)
    # Prepare the templates and emails
    email_template = "herald/herald_email.html"
    author = post.author.profile.name
    from_email = settings.DEFAULT_NOREPLY_EMAIL

    # Total number of recipients allowed per open connection,
    # Amazon SES has limit of 50
    batch_size = 40

    # Iterate through recipients and send emails in batches.
    for idx in range(0, len(emails), batch_size):
        # Get the next set of emails
        end = idx + batch_size
        rec_list = emails[idx:end]
        send_email(template_name=email_template, extra_context=context, name=author,
        from_email=from_email,recipient_list=rec_list, mass=True)

    return


@task
def mailing_list(uid, extra_context={}):
    """
    Generate notification for mailing list users.
    """
    from django.conf import settings
    from biostar.forum.models import Post
    from biostar.accounts.models import User, Profile

    # Get the post and users that have this enabled.
    post = Post.objects.filter(uid=uid).first()
    users = User.objects.filter(profile__digest_prefs=Profile.ALL_MESSAGES)

    emails = [user.email for user in users]

    # Update template context with post
    extra_context.update(dict(post=post))

    # Prepare the templates and emails
    email_template = "messages/mailing_list.html"
    author = post.author.profile.name
    from_email = settings.DEFAULT_NOREPLY_EMAIL
    if emails:
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
    subs = Subscription.objects.filter(id__in=sub_ids)

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
    email_subs = subs.filter(type=Subscription.EMAIL_MESSAGE)
    # Exclude mailing list users to avoid duplicate emails.
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
