import logging
from django.conf import settings
from biostar.message import models, auth
from django.contrib.auth import get_user_model

User = get_user_model()

logger = logging.getLogger("engine")

HAS_UWSGI = False

COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    @spool(pass_arguments=True)
    def async_create_sub_messages(subs, author, content, post_type):
        """
        Create messages to users subscribed to a post.
        """

        return


    @spool(pass_arguments=True)
    def async_notify_mentions(users, root, author, content):
        """
        Create local messages when users get mentioned in a post
        """

        return

except ModuleNotFoundError as exc:
    pass


def notify_mentions(users, root, author, content):
    title = root.title
    body = f"""
            Hello,

            You have been mentioned in a post by {author.profile.name}.

            The root post is :{title}.

            Here is where you are mentioned :

            {content}
            """
    subject = f"Mentioned in a post: {title}"

    message = auth.create_messages(body=body, subject=subject, recipient_list=users,
                                   mtype=models.Message.LOCAL_MESSAGE, sender=author,
                                   source=models.Message.MENTIONED)

    return message


def create_discussion_messages():
    return


def create_sub_messages(subs, root, author, content):
    "Create subscription messages"

    users_id_list = subs.values_list("user", flat=True).distinct()
    title = root.title

    subbed_users = User.objects.filter(id__in=users_id_list)

    body = f"""

        Hello,\n

        There is an addition by {author.profile.name} to a post you are subscribed to.\n

        Post: {title}\n

        Addition: {content}\n

        """

    subject = f"Subscription to : {title}"

    message = auth.create_messages(body=body, subject=subject, recipient_list=subbed_users,
                                   mtype=models.Message.LOCAL_MESSAGE, sender=author)

    return message
