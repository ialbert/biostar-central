from datetime import timedelta
import logging
from django.template import loader
from django.core.management.base import BaseCommand
from taggit.models import Tag
from biostar.forum.models import Post
from biostar.emailer.tasks import send_email
from biostar.accounts import util, models

logger = logging.getLogger('biostar')


def send_digests(days=1, subject=""):
    '''
    Send digest emails to users
    '''

    mapper = {1: models.Profile.DAILY_DIGEST,
              7: models.Profile.WEEKLY_DIGEST,
              30: models.Profile.MONTHLY_DIGEST}

    # Get posts made within the given time range.
    trange = util.now() - timedelta(days=days)

    posts = Post.objects.filter(lastedit_date__gt=trange, is_toplevel=True).order_by('-lastedit_date')

    if not posts:
        logger.info(f'No new posts found in the last {days} days.')
        return

    # Get users with the appropriate digest preference.
    pref = mapper.get(days, models.Profile.DAILY_DIGEST)
    users = models.User.objects.filter(profile__digest_prefs=pref)

    # Prepare the email
    context = dict(subject=subject, posts=posts)
    # Queue and send digest emails.
    emails = users.values_list('email', flat=True)
    send_email(template_name="messages/digest.html", extra_context=context, recipient_list=emails)

    return


class Command(BaseCommand):
    help = 'Send user digests to users.'

    def add_arguments(self, parser):

        parser.add_argument('--daily', dest='daily', action='store_true', help='Send daily digests.')
        parser.add_argument('--weekly', dest='weekly', action='store_true', help='Send weekly digests.')
        parser.add_argument('--monthly', dest='monthly', action='store_true', help='Send monthly digests.')

    def handle(self, *args, **options):
        monthly = options['monthly']
        weekly = options['weekly']
        daily = options['daily']

        if daily:
            send_digests(days=1, subject="Daily digest")
        elif weekly:
            send_digests(days=7, subject="Weekly digest")
        elif monthly:
            send_digests(days=30, subject="Monthly digest")
