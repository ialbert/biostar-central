from __future__ import  unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand
from django.conf import settings
from django.core.mail import send_mail
from biostar.const import now
from datetime import datetime, timedelta
logger = logging.getLogger(__name__)
from optparse import make_option
from biostar.apps.util import html
from django.contrib.sites.models import Site

def daily_digest(start=None, send=False):
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User

    site = site = Site.objects.get_current()

    start = start or (now() - timedelta(hours=24))

    # Select the 24 hour period since the last email.
    end = start + timedelta(hours=24)

    created = Post.objects.filter(status=Post.OPEN, type__in=Post.TOP_LEVEL, creation_date__gt=start, creation_date__lt=end).select_related('author')

    updated = Post.objects.filter(status=Post.OPEN, lastedit_date__gt=start, lastedit_date__lt=end)
    updated = updated.exclude(creation_date__gt=start, creation_date__lt=end, type__in=Post.TOP_LEVEL).select_related("author")

    blogs = Post.objects.filter(status=Post.OPEN, type=Post.BLOG, creation_date__gt=start, creation_date__lt=end).select_related('author')

    post_count = Post.objects.filter(status=Post.OPEN).count()

    params = dict(
        site = site,
        created=created,
        updated=updated,
        post_count=post_count
    )
    html_body = html.render("explorer/daily_digest.html", **params)

    print (html_body)

    # params = dict(subject=subject, from_email=from_email, recipient_list=recipient_list)
    #send_mail(subject=subject, message=message, from_email=from_email, recipient_list=recipient_list)


class Command(BaseCommand):
    help = 'send digest emails'

    option_list = BaseCommand.option_list + (

        make_option('--weekly', dest='weekly', action='store_true', default=False,
                    help='send the email to the weekly email subscribers'),

        make_option('--send', dest='send', action='store_true', default=False,
                    help='send the email'),

        make_option('-f', dest='file', default=None,
                    help='read email content from file (should be a template)'),

    )

    def handle(self, *args, **options):
        from_email = settings.DEFAULT_FROM_EMAIL

        if options['weekly']:
            pass
        else:
            daily_digest()




