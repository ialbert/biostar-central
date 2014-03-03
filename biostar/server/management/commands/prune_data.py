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
    from biostar.apps.posts.models import PostView
    from biostar.apps.messages.models import Message

    # Reduce post views.
    past = now() - timedelta(days=days)
    query = PostView.objects.filter(date__lt=past)
    msg = "Deleting %s PostViews" % query.count()
    logger.info(msg)
    query.delete()

    # Reduce messages.
    since = now() - timedelta(weeks=weeks)
    query = Message.objects.filter(sent_at__lt=since)
    msg = "Deleting %s messages" % query.count()
    logger.info(msg)
    query.delete()

if __name__ == '__main__':
    #generate_sitemap()
    pass