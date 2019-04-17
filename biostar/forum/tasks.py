import logging
from biostar.accounts.models import User
from biostar.forum.models import Award, Badge, Post
from biostar.forum.awards import ALL_AWARDS
from biostar.forum import util
from biostar.message.models import Message
#from biostar.emailer.auth import notify
from django.core.mail import send_mass_mail

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

    # 3600 secs in an hour X 24 hours in a day.
    secs = days * 3600 * 24

    return secs


try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    def send_post_mail(pid):
        """Send mail for users subscribed to specific post """

        return

    @timer(secs=days_to_secs(7))
    def send_weekly_digest():
        """Send weekly digest to users """

        # Get top level posts a week old

        return

    @timer(secs=days_to_secs())
    def send_daily_digest():
        """Send daily digest to users """

        # Get top level post created today


        return

    @timer(secs=days_to_secs(30))
    def send_monthly_digest():
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
                    context = '<a href="%s">%s</a>' % (target.get_absolute_url(), target.title)
                    award = Award.objects.create(user=user, badge=badge, date=date, post=target,
                                                 context=context)
                else:
                    context = ""
                    award = Award.objects.create(user=user, badge=badge, date=date, context=context)
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
