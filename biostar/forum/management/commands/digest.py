from datetime import timedelta
from django.template import loader
from django.core.management.base import BaseCommand
from taggit.models import Tag
from biostar.forum.models import Post
from biostar.emailer.tasks import send_email
from biostar.accounts import util, models


def send_digests(days=1, subject=""):
    '''
    Send digest emails to users
    '''
    prefs_map = {1: models.Profile.DAILY_DIGEST, 7: models.Profile.WEEKLY_DIGEST, 30: models.Profile.MONTHLY_DIGEST}
    # Send digests for posts within the last x days.
    delta = util.now() - timedelta(days=days)

    # Get users with the appropriate digest prefs.
    digest_prefs = prefs_map.get(days, models.Profile.DAILY_DIGEST)
    users = models.User.objects.filter(profile__digest_prefs=digest_prefs)

    # Fetch posts within the last x amount of days
    posts = Post.objects.filter(lastedit_date__gt=delta)

    email_template = loader.get_template("messages/digest.html")
    context = dict(subject=subject, posts=posts)

    # Queue and send digest emails.
    emails = users.values_list('email', flat=True)
    send_email(template_name=email_template, extra_context=context, recipient_list=emails)

    return


def send_watched_tags():
    '''
    Get watched tags from the last week and send email to users.
    '''

    # Get users with watched tags.
    users = models.User.objects.exclude(profile__watched_tags='')

    # Get posts made in the last week
    delta = util.now() - timedelta(days=7)

    posts = Post.objects.filter(lastedit_date__lte=delta)

    tags = Tag.objects.filter(post__lastedit_date__gt=delta, post__type__in=Post.TOP_LEVEL)
    tags = tags.values_list('name', flat=True)

    # Fetch users with the watched tags
    #TODO: this will probably need to be done in a for loop.
    users = users.filter(watched_tags__icontains=tags)

    print(users)

    return


def test():
    # Testing the cron tasks
    print("in main")

    print("TESTING")

    import os
    test = os.path.abspath(os.path.join(os.path.dirname(__file__), 'tmp'))
    print(test)

    open(test, 'w').write('TESTING this cron')

    return


class Command(BaseCommand):
    help = 'Send user digests to users.'

    def add_arguments(self, parser):

        parser.add_argument('--daily', dest='daily', action='store_true', help='Send daily digests.')
        parser.add_argument('--weekly', dest='weekly', action='store_true', help='Send weekly digests.')
        parser.add_argument('--monthly', dest='monthly', action='store_true', help='Send monthly digests.')
        parser.add_argument('--tags', dest='tags', action='store_true', help='Send watched tags digests to users.')

    def handle(self, *args, **options):
        tags = options['tags']
        monthly = options['monthly']
        weekly = options['weekly']
        daily = options['daily']

        # Index all un-indexed posts that have a root.
        if daily:
            send_digests(days=1, subject="Daily digest")
        elif weekly:
            send_digests(days=7, subject="Weekly digest")
        elif monthly:
            send_digests(days=30, subject="Monthly digest")
        elif tags:
            send_watched_tags()
