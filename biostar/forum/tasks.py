import logging
from django.conf import settings
logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @spool(pass_arguments=True)
    def create_messages(subs):
        """
        Create messages to users subscribed to a post.
        """
        return

    @spool(pass_arguments=True)
    def notify_mentions(users, post):
        """
        Create local messages when users get mentioned in a post
        """

        return



except ModuleNotFoundError as exc:
    pass
