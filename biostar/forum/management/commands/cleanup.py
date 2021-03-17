
import logging
from datetime import datetime, timedelta
from django.db.models import Count
from django.core.management.base import BaseCommand
from biostar.accounts.models import Message, User
from biostar.forum.util import now
from biostar.forum.models import PostView, Post

logger = logging.getLogger('engine')


MAX_MSG = 100


def prune_data(weeks=10, days=1, delall=False):

    # Remove spam and deleted posts
    past_days = now() - timedelta(days=days)
    weeks_since = now() - timedelta(weeks=weeks)

    #spam_posts = Post.objects.filter(spam=Post.SPAM)
    #logger.info(f"Deleting {spam_posts.count()} spam posts")
    #spam_posts.delete()

    # Remove post views.
    post_views = PostView.objects.filter(date__lt=past_days)
    logger.info(f"Deleting {post_views.count()} post views")
    post_views.delete()

    # Reduce overall messages.
    if delall:
        messages = Message.objects.all()
    else:
        messages = Message.objects.filter(sent_date__lt=weeks_since)

    logger.info(f"Deleting {messages.count()} messages")
    messages.delete()

    # Get rid of too many messages
    #users = User.objects.annotate(total=Count("message__recipient")).filter(total__gt=MAX_MSG)[:100]
    #for user in users:
    #    since = now() - timedelta(days=1)
    #    Message.objects.filter(user=user, sent_at__lt=since).delete()

    return


class Command(BaseCommand):
    help = """Delete the following: 
              - posts marked as spam
              - messages older then 10 weeks
              - PostView objects within the last day 
              - too many messages in a users inbox
           """

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', default=False, help="Delete all messages.")

    def handle(self, *args, **options):

        delete_all = options['all']

        prune_data(delall=delete_all)