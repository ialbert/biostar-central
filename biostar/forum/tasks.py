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


def send_message(template, context, sender, rec_list):
    from django.template import loader
    from .auth import create_local_messages
    # Render the template
    tmpl = loader.get_template(template_name=template)
    body = tmpl.render(context)

    # Send the local message
    create_local_messages(body=body, sender=sender, rec_list=rec_list)


def check_profile(request, user):
    from biostar.accounts import auth
    auth.check_user_profile(request=request, user=user)
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