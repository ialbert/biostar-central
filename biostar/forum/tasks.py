import logging
from biostar.accounts.models import User
from biostar.forum.models import Award, Badge, Post
from biostar.forum.awards import ALL_AWARDS

logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    def send_mail():
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

                if isinstance(target, Post):
                    context = '<a href="%s">%s</a>' % (target.get_absolute_url(), target.title)
                else:
                    context = ""

                date = user.profile.last_login
                award = Award.objects.create(user=user, badge=badge, date=date, context=context)
                logger.info("award %s created for %s" % (award.badge.name, user.email))
                print(target, award)
        return


    def check_user_profile(ip, user):
        return

except ModuleNotFoundError as exc:
    pass