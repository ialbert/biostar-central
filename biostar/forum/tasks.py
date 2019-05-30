import logging
from django.conf import settings
from biostar.accounts import tasks, models
from biostar.forum.models import Subscription
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


def send_email(template, context, subject, email_list, from_email):
    from biostar.emailer.auth import send
    send(template_name=template, email_list=email_list,
         extra_context=context, from_email=from_email,
         subject=subject, send=True)


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
    subs = Subscription.objects.filter(post=post.root).exclude(type=models.Profile.NO_MESSAGES)

    # Send local messages
    users = [sub.user for sub in subs if sub.user != author]
    tasks.create_messages.spool(template=local_template, extra_context=context, rec_list=users, sender=author)

    # Send emails to users that specified so
    subs = subs.filter(type=models.Profile.EMAIL_MESSAGE)
    emails = [sub.user.email for sub in subs if (sub.user != author and sub.type == models.Profile.EMAIL_MESSAGE)]

    from_email = settings.ADMIN_EMAIL

    send_email.spool(template=email_template, context=context, subject="Subscription",
                     email_list=emails, from_email=from_email)


if HAS_UWSGI:
    info_task = spool(info_task, pass_arguments=True)
    send_email = spool(send_email, pass_arguments=True)
    notify_followers = spool(notify_followers, pass_arguments=True)
    created_post = spool(created_post, pass_arguments=True)
else:
    info_task.spool = info_task
    notify_followers.spool = notify_followers
    send_email.spool = send_email
    created_post.spool = created_post
