import logging
from datetime import datetime
from django.utils.timezone import utc
from biostar.accounts.models import User, Profile
from biostar.forum.models import Award, Badge, Post
from biostar.forum.awards import ALL_AWARDS
from biostar.forum import util
from biostar.message.models import Message
from django.template import loader

from django.conf import settings
from biostar.emailer.sender import EmailTemplate
from django.core.mail import send_mail

logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1


def send_award_messages(award):
    """ Send local message to user when award is won"""

    user = award.user

    award_template = "message/award_creation.html"
    # Generate the message from the template.
    content = util.render(name=award_template, award=award, user=user)

    subject = "Congratulations: you won %s" % award.badge.name

    # Create the message body.
    Message.objects.create(user=user, subject=subject, body=content,)

    return


def days_to_secs(days=1):
    """
    Convert days to seconds
    """
    # 3600 secs in an hour X 24 hours in a day.
    secs = days * 3600 * 24
    return secs


def send_digest_emails(template_name, extra_context, digest_prefs=Profile.WEEKLY_DIGEST):

    # Get users that opted for digest
    users = User.objects.filter(profile__digest_prefs=digest_prefs,
                                profile__message_prefs=Profile.DIGEST_MESSAGES)
    # debug only
    to_email = users.filter(email="natay.aberra@gmail.com").first()
    to_email = to_email.values_list("email", flat=True)
    #to_email = users.values_list("email", flat=True)

    # Add site info to context
    context = dict(domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL, port=settings.HTTP_PORT, name=settings.SITE_NAME)
    context.update(extra_context)

    # Load the template and prepare email
    template = loader.get_template(template_name)
    email = template.render(context=context)
    subject = context.get("subject", "Digest")

    print(to_email, email)
    1 / 0
    # Get the from email
    from_email = User.objects.filter(is_superuser=True).first()
    # Render the template

    # Send email.
    send_mail(subject=subject, from_email=from_email, recipient_list=rec_list, html_message=email, message=email)
    return


try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    def send_post_mail(pid):
        """Send email for users subscribed to specific post """

        return

    #@timer(secs=1)
    @timer(secs=days_to_secs(days=1))
    def send_daily_digest(args):
        """Send daily digest to users """
        today = datetime.utcnow().replace(tzinfo=utc)
        posts = Post.objects.filter(type__in=Post.TOP_LEVEL, creation_date=today)

        # Prepare context to populate template
        subject = f"Daily Digest for :{today.date()}"
        msg = f"Hello, here are a digest of posts from today, {today.date()}. "
        context = dict(posts=posts, msg=msg, subject=subject)

        #send_digest_emails(digest_prefs=Profile.DAILY_DIGEST,
        #                   extra_context=context, template_name="messages/digest.html")

        return

    @timer(secs=days_to_secs(days=7))
    def send_weekly_digest(args):
        """Send weekly digest to users """

        today = datetime.utcnow().replace(tzinfo=utc)
        posts = Post.objects.filter(type__in=Post.TOP_LEVEL, creation_date=today)

        # Load template with message
        subject = f"weekly Digest for :{today.date()}"
        msg = f"Hello, here are a digest of posts from today, {today.date()}. "
        context = dict(posts=posts, msg=msg, subject=subject)

        #send_digest_emails(template_name="messages/digest.html", digest_prefs=Profile.WEEKLY_DIGEST,
        #                   extra_context=context)

        return


    @timer(secs=days_to_secs(days=30))
    def send_monthly_digest(args):
        """Send monthly (30 days) digest to users """

        # Get top level posts one month old.

        return


    @spool(pass_arguments=True)
    def create_user_awards(user_id):

        user = User.objects.filter(id=user_id).first()
        awards = dict()
        for award in Award.objects.filter(user=user).select_related('badge'):
            awards.setdefault(award.badge.name, []).append(award)

        get_award_count = lambda name: len(awards[name]) if name in awards else 0

        for obj in ALL_AWARDS:

            # How many times has been awarded to this user.
            seen_count = get_award_count(obj.name)

            # How many times should it been awarded
            valid_targets = obj.validate(user)

            # Keep that targets that have not been awarded
            valid_targets = valid_targets[seen_count:]

            # Some limit on awards
            valid_targets = valid_targets[:100]

            # Award the targets
            for target in valid_targets:
                # Update the badge counts.
                badge = Badge.objects.filter(name=obj.name)
                badge.count += 1
                badge.save()
                date = user.profile.last_login

                if isinstance(target, Post):
                    award = Award.objects.create(user=user, badge=badge, date=date, post=target)
                else:
                    award = Award.objects.create(user=user, badge=badge, date=date)

                send_award_messages(award=award)
                logger.info("award %s created for %s" % (award.badge.name, user.email))

        return

    def check_user_profile(ip, user):
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


except ModuleNotFoundError as exc:
    pass
