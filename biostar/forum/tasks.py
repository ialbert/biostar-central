import logging

from django.conf import settings
from django.db.models import F, Q
from django.shortcuts import reverse
from biostar.accounts.tasks import create_messages, detect_location
from biostar.emailer.tasks import send_email
from biostar.utils.decorators import spool

logger = logging.getLogger("biostar")


@spool(pass_arguments=True)
def info_task(*args, **kwargs):
    logger.info(f"info_task called with args={args} and kwargs={kwargs}")


@spool(pass_arguments=True)
def created_post(pid):
    logger.info(f"Created post={pid}")


def create_award(targets, user, award):
    from biostar.forum.models import Award, Post, Badge

    for target in targets:
        date = user.profile.last_login
        post = target if isinstance(target, Post) else None
        badge = Badge.objects.filter(name=award.name).first()
        Award.objects.create(user=user, badge=badge, date=date, post=post)

        logger.debug("award %s created for %s" % (badge.name, user.email))

    return


@spool(pass_arguments=True)
def create_user_awards(user_id):
    from biostar.accounts.models import Profile, User
    from biostar.forum.models import Award, Post
    from biostar.forum.awards import ALL_AWARDS

    user = User.objects.filter(id=user_id).first()
    if (user.profile.state == Profile.NEW) and (user.profile.score > 10):
        user.profile.state = Profile.TRUSTED
        user.save()
    # The awards the user has won at this point
    awards = dict()
    for award in Award.objects.filter(user=user).select_related('badge'):
        awards.setdefault(award.badge.name, []).append(award)

    for award in ALL_AWARDS:

        # How many times has this award has already been given
        seen = len(awards[award.name]) if award.name in awards else 0

        # How many times the user earned this award
        valid_targets = award.validate(user)

        # Keep targets have not been awarded
        valid_targets = valid_targets[seen:]

        # Create an award for each target
        create_award(targets=valid_targets, user=user, award=award)


@spool(pass_arguments=True)
def notify_followers(post, author):
    """
    Send subscribed users, excluding author, a message/email.
    """
    from biostar.forum.models import Subscription

    # Template used to send local messages
    local_template = "messages/subscription_message.md"

    # Template used to send emails with
    email_template = "messages/subscription_email.html"

    # Everyone subscribed gets a local message.
    subs = Subscription.objects.filter(post=post.root).exclude(Q(type=Subscription.NO_MESSAGES) | Q(user=author))

    # Does the does not have subscriptions.
    if not subs:
        return

    # Select users that should be notified.
    users = [sub.user for sub in subs]

    # Additional context for the message.
    extra_context = dict(post=post)

    # Every use gets local messages if subscribed in any way.
    create_messages(template=local_template, extra_context=extra_context, rec_list=users, sender=author)

    # Select users with email subscriptions.
    email_subs = subs.filter(type=Subscription.EMAIL_MESSAGE)

    # No email subscriptions
    if not email_subs:
        return

    emails = [sub.user.email for sub in subs]

    send_email(template_name=email_template, extra_context=extra_context, subject="Subscription",
               email_list=emails)
