import logging
import re
from django.conf import settings
from biostar.message import models, auth
from biostar.accounts.models import Profile
from biostar.forum.models import Subscription
from django.contrib.auth import get_user_model

User = get_user_model()

logger = logging.getLogger("engine")

HAS_UWSGI = False

COUNTER = 1


def parse_mentioned_users(content):

    # Any word preceded by a @ is considered a user handler.
    handler_pattern = "\@[^\s]+"
    # Drop leading @
    users_list = set(x[1:] for x in re.findall(handler_pattern, content))

    return User.objects.filter(username__in=users_list)


try:

    from uwsgidecorators import *

    HAS_UWSGI = True

    @spool(pass_arguments=True)
    def async_create_messages(subject, author, body, rec_list, source=models.Message.MENTIONED):
        """
        Create messages to users in recipient list
        """
        # Assign a task to a a worker
        auth.create_local_messages(body=body, subject=subject, rec_list=rec_list, sender=author, source=source)


except (ModuleNotFoundError, NameError) as exc:
    HAS_UWSGI = False
    # Bail out and record error
    logger.error(exc)
    pass


def parse_mention_msg(post, author):
    title = post.title
    rec_list = parse_mentioned_users(content=post.content)
    defalut_body = f"""
            Hello, You have been mentioned in a post by {author.profile.name}.
            The root post is :{title}.
            Here is where you are mentioned :
            {post.content}
            """
    body = settings.MENTIONED_TEMPLATE or defalut_body
    subject = f"Mentioned in a post: {title}"

    return body, subject, rec_list


def parse_subs_msg(post, author):
    title = post.title
    # Load template if its available
    default_body = f"""
          Hello,\n
          There is an addition by {author.profile.name} to a post you are subscribed to.\n
          Post: {title}\n
          Addition: {post.content}\n
          """
    root = post if post.is_toplevel else post.root
    rec_list = Subscription.objects.filter(post=root).user_set()

    body = settings.SUBS_TEMPLATE or default_body
    subject = f"Subscription to : {title}"
    return body, subject, rec_list


def create_message(post, author):

    ment_subject, ment_body, ment_users = parse_mention_msg(post=post, author=author)
    sub_subject, sub_body, sub_users = parse_subs_msg(post=post,  author=author)

    # Create asynchronously when uwsgi is available
    if HAS_UWSGI:
        # Assign a worker to send subscription messages
        async_create_messages(subject=sub_subject, body=sub_body, rec_list=sub_users)

        # Assign a worker to send mentioned users
        async_create_messages(source=models.Message.MENTIONED,
                              subject=ment_subject, body=ment_body, rec_list=ment_users)

        return

    # Can run synchrony only when debugging
    if settings.DEBUG:
        # Send mentioned messages
        auth.create_local_messages(body=ment_body, subject=ment_subject, rec_list=sub_users,
                                   sender=author, source=models.Message.MENTIONED)
        # Send subscription messages
        auth.create_local_messages(body=sub_body, subject=sub_subject, rec_list=ment_users,
                                   sender=author)

        return

