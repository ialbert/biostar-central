from __future__ import absolute_import
from django.conf import settings

from .celery import app

import logging

logger = logging.getLogger(__name__)


def init_awards():
    "Initializes the badges"
    from biostar.apps.badges.models import Badge
    from biostar.apps.badges.award_defs import ALL_AWARDS

    for obj in ALL_AWARDS:
        badge, created = Badge.objects.get_or_create(name=obj.name)

        # Badge descriptions may change.
        if badge.desc != obj.desc:
            badge.desc = obj.desc
            badge.icon = obj.icon
            badge.type = obj.type
            badge.save()

        if created:
            logger.info("initializing badge %s" % badge)


@app.task
# Tries to award a badge to the user
def create_post_award(user):
    pass


@app.task
# Tries to award a badge to the user
def create_user_award(user):
    from biostar.apps.users.models import User
    from biostar.apps.badges.models import Badge, Award
    from biostar.apps.badges.award_defs import ALL_AWARDS

    # Update user types.
    if (user.type == User.NEW_USER) and (user.score > 10):
        user.type = User.TRUSTED
        user.save()

    # Debug only
    #Award.objects.all().delete()

    # The awards the user has won at this point
    awards = dict()
    for award in Award.objects.filter(user=user).select_related('badge'):
        awards.setdefault(award.badge.name, []).append(award)

    # Function to detect if the award already exists
    get_award_count = lambda name: len(awards[name]) if name in awards else 0

    for obj in ALL_AWARDS:

        # How many times has been awarded
        seen_count = get_award_count(obj.name)

        # How many times should it been awarded
        valid_count = obj.validate(user)

        if valid_count > seen_count:

            # Update the badge counts.
            badge = Badge.objects.get(name=obj.name)
            badge.count += (valid_count - seen_count)
            badge.save()

            # Create the corresponding awards as many times as needed.
            for x in range(valid_count - seen_count):
                award = Award.objects.create(user=user, badge=badge)
                logger.info("award %s created for %s" % (award.badge.name, user.email))

            # Grant one type of award at a time.
            break