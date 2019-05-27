import logging
from urllib.request import urlopen
import json
from biostar.accounts.models import User, Profile
from biostar.accounts import auth as accounts_auth

from .auth import create_local_messages

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
    def async_create_messages(subject, sender, body, rec_list, uid=None):
        """
        Create messages to users in recipient list
        """
        # Assign a task to a a worker
        create_local_messages(body=body, subject=subject, rec_list=rec_list, sender=sender, uid=uid)

    @spool(pass_arguments=True)
    def async_check_profile(request, user_id):
        user = User.objects.filter(id=user_id)
        accounts_auth.check_user_profile(request=request, user=user)
        logger.info(f"Checked user profile user={user}")

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


def send_message(subject, body, rec_list, sender, uid=None):
    # Create asynchronously when uwsgi is available
    if HAS_UWSGI:
        # Assign a worker to send mentioned users
        async_create_messages(sender=sender, subject=subject, body=body, rec_list=rec_list, uid=uid)
        return
    # Can run synchrony only when debugging
    if settings.DEBUG:
        # Send subscription messages
        auth.create_local_messages(body=body, sender=sender, subject=subject, rec_list=rec_list, uid=uid)
    return