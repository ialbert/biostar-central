from __future__ import absolute_import
from django.conf import settings

from .celery import app

import logging


from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)

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
def check_user_profile(ip, user):
    import urllib2, json
    logger.info("profile check from %s on %s" % (ip, user))
    if not user.profile.location:
        try:
            url = "http://api.hostip.info/get_json.php?ip=%s" % ip
            logger.info("%s, %s, %s" % (ip, user, url))
            f = urllib2.urlopen(url, timeout=3)
            data = json.loads(f.read())
            f.close()
            location = data.get('country_name', '').title()
            if "unknown" not in location.lower():
                user.profile.location = location
                user.profile.save()
        except Exception, exc:
            logger.error(exc)

@app.task
# Tries to award a badge to the user
def create_user_award(user):
    from biostar.apps.users.models import User
    from biostar.apps.posts.models import Post
    from biostar.apps.badges.models import Badge, Award
    from biostar.apps.badges.award_defs import ALL_AWARDS

    logger.info("award check for %s" % user)

    # Update user status.
    if (user.status == User.NEW_USER) and (user.score > 10):
        user.status = User.TRUSTED
        user.save()

    # Debug only
    #Award.objects.all().delete()

    # The awards the user has won at this point
    awards = dict()
    for award in Award.objects.filter(user=user).select_related('badge'):
        awards.setdefault(award.badge.name, []).append(award)

    # Shorcut function to get the award count
    get_award_count = lambda name: len(awards[name]) if name in awards else 0

    for obj in ALL_AWARDS:

        # How many times has been awarded
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
            badge = Badge.objects.get(name=obj.name)
            badge.count += 1
            badge.save()

            if isinstance(target, Post):
                context = '<a href="%s">%s</a>' % (target.get_absolute_url(), target.title)
            else:
                context = ""

            date = user.profile.last_login
            award = Award.objects.create(user=user, badge=badge, date=date, context=context)
            logger.info("award %s created for %s" % (award.badge.name, user.email))

