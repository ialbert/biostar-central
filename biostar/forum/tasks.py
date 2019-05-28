import logging
from urllib.request import urlopen
import json
from biostar.accounts.models import User, Profile
from biostar.accounts import auth as accounts_auth

from biostar.emailer.auth import notify
from biostar.message.auth import create_local_messages


logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1


def days_to_secs(days=1):
    """
    Convert days to seconds
    """
    # 3600 secs in an hour X 24 hours in a day.
    secs = days * 3600 * 24
    return secs


try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @spool(pass_arguments=True)
    def async_check_profile(request, user_id):
        user = User.objects.filter(id=user_id)
        accounts_auth.check_user_profile(request=request, user=user)
        logger.info(f"Checked user profile user={user}")

    @spool(pass_arguments=True)
    def async_create_messages(subject, sender, body, rec_list, parent=None, uid=None):
        """
        Create messages to users in recipient list
        """
        # Assign a task to a a worker
        create_local_messages(body=body, subject=subject, rec_list=rec_list, sender=sender, parent=parent, uid=uid)


    @spool(pass_arguments=True)
    def async_send_email(template, emails, context, from_email, subject, send=False):
        notify(template_name=template, email_list=emails,
               extra_context=context, from_email=from_email,
               subject=subject, send=True)

        # template = "default_messages/subscription_msg.html"
        #
        # context = dict(post=post)
        # tmpl = loader.get_template(template_name=template)
        # body = tmpl.render(context)
        #
        # subs = Subscription.objects.filter(post=post.root)
        # subs = subs.exclude(type=Profile.NO_MESSAGES)
        # id_list = subs.values_list("user", flat=True).exclude(id=author.pk).distinct()

        return

    @spool(pass_arguments=True)
    def created_post(pid):
        logger.info(f"Created post={pid}")

    @spool(pass_arguments=True)
    def edited_post(pid):
        logger.info(f"Edited post={pid}")

    @spool(pass_arguments=True)
    def added_sub(sid):
        logger.info(f"Created sub with pk={sid}")

    @spool(pass_arguments=True)
    def moderated_post(pid):
        logger.info(f"Post has been moderated pid={pid}")

    @spool(pass_arguments=True)
    def triggered_vote(pid, vtype):
        logger.info(f"Created Vote for post={pid} with type={vtype}")


except (ModuleNotFoundError, NameError) as exc:
    HAS_UWSGI = False
    logger.error(exc)
    pass


def send_message(subject, body, rec_list, sender, parent=None, uid=None):
    # Create asynchronously when uwsgi is available
    if HAS_UWSGI:
        # Assign a worker to send mentioned users
        async_create_messages(sender=sender, subject=subject, body=body,
                              rec_list=rec_list, parent=parent, uid=uid)
        return

    # Send messages synchronously
    create_local_messages(body=body, sender=sender, subject=subject, rec_list=rec_list,
                          parent=parent, uid=uid)
    return


def send_email():
    return





