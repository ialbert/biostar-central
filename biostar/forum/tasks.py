import logging

from django.db.models import Q

from biostar.accounts.tasks import create_messages
from biostar.emailer.tasks import send_email
from biostar.utils.decorators import spool

logger = logging.getLogger("biostar")


@spool(pass_arguments=True)
def info_task(*args, **kwargs):
    logger.info(f"info_task called with args={args} and kwargs={kwargs}")


@spool(pass_arguments=True)
def created_post(pid):
    logger.info(f"Created post={pid}")


@spool(pass_arguments=True)
def create_user_awards(user_id):
    from biostar.accounts.models import User
    from biostar.forum.models import Award, Badge, Post
    from biostar.forum.awards import ALL_AWARDS

    user = User.objects.filter(id=user_id).first()

    # debugging
    #Award.objects.all().delete()

    for award in ALL_AWARDS:
        # Valid award targets the user has earned
        targets = award.validate(user)

        for target in targets:
            date = user.profile.last_login
            post = target if isinstance(target, Post) else None
            badge = Badge.objects.filter(name=award.name).first()

            # Do not award a post multiple times.
            already_awarded = Award.objects.filter(user=user, badge=badge, post=post).exists()
            if post and already_awarded:
                continue

            # Create an award for each target.
            Award.objects.create(user=user, badge=badge, date=date, post=post)

            logger.debug("award %s created for %s" % (badge.name, user.email))


@spool(pass_arguments=True)
def notify_followers(subs, author, extra_context={}):
    """
    Generate notification to users subscribed to a post, excluding author, a message/email.
    """
    from biostar.forum.models import Subscription

    # Template used to send local messages
    local_template = "messages/subscription_message.md"

    # Template used to send emails with
    email_template = "messages/subscription_email.html"

    # Does the does not have subscriptions.
    if not subs:
        return

    # Select users that should be notified.
    users = [sub.user for sub in subs]

    # Every subscribed user gets local messages with any subscription type.
    create_messages(template=local_template, extra_context=extra_context, rec_list=users, sender=author)

    # Select users with email subscriptions.
    email_subs = subs.filter(type=Subscription.EMAIL_MESSAGE)

    # No email subscriptions
    if not email_subs:
        return

    recipient_list = [sub.user.email for sub in subs]

    send_email(template_name=email_template, extra_context=extra_context, recipient_list=recipient_list)
