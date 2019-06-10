import logging

from django.conf import settings
from django.db.models import F
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
    template = "messages/awards_created.md"

    for target in targets:
        date = user.profile.last_login
        post = target if isinstance(target, Post) else None
        badge = Badge.objects.filter(name=award.name).first()
        awarded = Award.objects.create(user=user, badge=badge, date=date, post=post)

        badge_url = reverse('badge_view', kwargs=dict(uid=badge.uid))
        context = dict(badge_url=badge_url, award=awarded, post=post)

        create_messages(template=template, extra_context=context, rec_list=[user])

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
    from biostar.accounts.models import Profile
    from biostar.forum.models import Subscription

    # Template used to send local messages
    local_template = "messages/subscription_message.md"
    # Template used to send emails with
    email_template = "messages/subscription_email.html"
    context = dict(post=post, post_url=post.get_absolute_url())

    # Everyone subscribed gets a local message.
    subs = Subscription.objects.filter(post=post.root).exclude(type=Profile.NO_MESSAGES)

    # Send local messages
    users = set(sub.user for sub in subs if sub.user != author)

    create_messages(template=local_template, extra_context=context, rec_list=users, sender=author)

    # Send emails to users that specified so
    subs = subs.filter(type=Profile.EMAIL_MESSAGE)
    emails = [sub.user.email for sub in subs if (sub.user != author and sub.type == Profile.EMAIL_MESSAGE)]

    from_email = settings.ADMIN_EMAIL

    send_email(template_name=email_template, extra_context=context, subject="Subscription",
               email_list=emails, from_email=from_email, send=True,)