import logging
from urllib.request import urlopen
import json


logger = logging.getLogger("engine")

HAS_UWSGI = False

try:
    from uwsgidecorators import *
    HAS_UWSGI = True

except (ModuleNotFoundError, NameError) as exc:
    HAS_UWSGI = False
    logger.error(exc)
    pass


def info_task(*args, **kwargs):
   logger.info(f"info_task called with {args} and {kwargs}")


def created_post(pid):
    logger.info(f"Created post={pid}")


def send_message(template, context, sender, subs=[]):
    from biostar.accounts import auth, models

    # Exclude current author of the post from receiving a message.
    user_ids = subs.values("user").exclude(user=sender).distinct()

    users = models.User.objects.filter(id__in=user_ids)
    # Send local messages
    auth.create_messages(template=template, sender=sender, rec_list=users, context=context)

    logger.debug(f"Sent to subscription message to {len(users)} users.")


def check_profile(request, user):
    from biostar.accounts.auth import check_user_profile
    check_user_profile(request=request, user=user)
    logger.info(f"Checked user profile user={user}")


def send_email(template, context, subject, email_list, from_email):
    from biostar.emailer.auth import notify
    notify(template_name=template, email_list=email_list,
           extra_context=context, from_email=from_email,
           subject=subject, send=True)


if HAS_UWSGI:
    info_task = spool(info_task, pass_arguments=True)
    send_message = spool(send_message, pass_arguments=True)
    check_profile = spool(send_message, pass_arguments=True)
    send_email = spool(send_email, pass_arguments=True)
    created_post = spool(created_post, pass_arguments=True)
else:
    info_task.spool = info_task
    send_message.spool = send_message
    check_profile.spool = check_profile
    send_email.spool = send_email
    created_post.spool = created_post