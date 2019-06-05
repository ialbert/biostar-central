import logging
from django.conf import settings
from django.db.models import F
from biostar.accounts.models import Profile, User
from biostar.accounts.tasks import detect_location, create_messages
from biostar.emailer.tasks import send_email
from biostar.forum.models import Subscription, Post, Award
logger = logging.getLogger("biostar")

try:
    from uwsgidecorators import *
    HAS_UWSGI = True
except (ModuleNotFoundError, NameError) as exc:
    HAS_UWSGI = False
    logger.warning(exc)
    pass


def info_task(*args, **kwargs):
   logger.info(f"info_task called with {args} and {kwargs}")


def created_post(pid):
    logger.info(f"Created post={pid}")


def create_user_awards(user):

    if (user.profile.state == Profile.NEW) and (user.score > 10):
        user.profile.state = Profile.TRUSTED
        user.save()

    # The awards the user has won at this point
    awards = dict()
    for award in Award.objects.filter(user=user).select_related('badge'):
        awards.setdefault(award.badge.name, []).append(award)

    print(awards)
    return


def create_subscription(root, user):

    # Create user subscription to post.
    sub, created = Subscription.objects.get_or_create(post=root, user=user)
    if created:
        # Increase subscription count of the root.
        Post.objects.filter(pk=root.pk).update(subs_count=F('subs_count') + 1)
        logger.debug(f"Created a subscription for user:{user} to root:{root.title}")

    return


def notify_followers(post, author):
    """
    Send subscribed users, excluding author, a message/email.
    """

    # Template used to send local messages
    local_template = "messages/subscription_message.html"
    # Template used to send emails with
    email_template = "messages/subscription_email.html"
    context = dict(post=post)

    # Everyone subscribed gets a local message.
    subs = Subscription.objects.filter(post=post.root).exclude(type=Profile.NO_MESSAGES)

    # Send local messages
    users = set(sub.user for sub in subs if sub.user != author)
    create_messages(template=local_template, extra_context=context, rec_list=users, sender=author)

    # Send emails to users that specified so
    subs = subs.filter(type=Profile.EMAIL_MESSAGE)
    emails = [sub.user.email for sub in subs if (sub.user != author and sub.type == Profile.EMAIL_MESSAGE)]

    from_email = settings.ADMIN_EMAIL

    send_email(template_name=email_template, extra_context=context, subject="Subscription",
               email_list=emails, from_email=from_email, send=True,)


if HAS_UWSGI:
    info_task = spool(info_task, pass_arguments=True)
    send_email = spool(send_email, pass_arguments=True)
    detect_location.spool = spool(detect_location, pass_arguments=True)
    create_messages.spool = spool(create_messages, pass_arguments=True)
    notify_followers = spool(notify_followers, pass_arguments=True)
    create_subscription = spool(create_subscription, pass_arguments=True)
    created_post = spool(created_post, pass_arguments=True)
    create_user_awards = spool(create_user_awards, pass_arguments=True)
else:
    info_task.spool = info_task
    create_messages.spool = create_messages
    detect_location.spool = detect_location
    notify_followers.spool = notify_followers
    create_subscription.spool = create_subscription
    send_email.spool = send_email
    created_post.spool = created_post
    create_user_awards.spool = create_user_awards
