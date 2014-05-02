"""
Prunes the views. Reduces
"""

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.utils.timezone import utc
from datetime import datetime, timedelta

import logging


def now():
    return datetime.utcnow().replace(tzinfo=utc)


logger = logging.getLogger("command")


class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def handle(self, *args, **options):
        main()


def main(days=1, weeks=10):
    from biostar.apps.posts.models import PostView, ReplyToken
    from biostar.apps.messages.models import Message
    from biostar.apps.users.models import User
    from django.db.models import Count

    # Reduce post views.
    past = now() - timedelta(days=days)
    query = PostView.objects.filter(date__lt=past)
    msg = "deleting %s views" % query.count()
    logger.info(msg)
    query.delete()

    # Reduce messages.
    since = now() - timedelta(weeks=weeks)
    query = Message.objects.filter(sent_at__lt=since)
    msg = "deleting %s messages" % query.count()
    logger.info(msg)
    query.delete()

    # Remove old reply tokens
    query = ReplyToken.objects.filter(date__lt=since)
    msg = "deleting %s tokens" % query.count()
    logger.info(msg)
    query.delete()

    # Get rid of too many messages
    MAX_MSG = 100
    users = User.objects.annotate(total=Count("recipients")).filter(total__gt=MAX_MSG)[:100]
    for user in users:
        since = now() - timedelta(days=1)
        Message.objects.filter(user=user, sent_at__lt=since).delete()


if __name__ == '__main__':
    #generate_sitemap()
    pass