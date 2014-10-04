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
from django.db.models import Count

def daily_digest(days, text_tmpl, html_tmpl, send, limit=10, verbosity=1):
    from_email = settings.DEFAULT_FROM_EMAIL
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User

    site = site = Site.objects.get_current()

    start = (now() - timedelta(days=days))

    # Posts created since the start date.
    top_posts = Post.objects.filter(status=Post.OPEN, type__in=Post.TOP_LEVEL,
                                  creation_date__gt=start).select_related('author')


    top_posts = top_posts.order_by('-view_count')[:limit]

    # Updated post created before the start date.
    upd_posts = Post.objects.filter(status=Post.OPEN,
                                  lastedit_date__gt=start)
    upd_posts = upd_posts.exclude(creation_date__gt=start, type__in=Post.TOP_LEVEL).select_related("author")
    upd_posts = upd_posts.order_by('-vote_count')[:limit]

    # Blog posts created since the start date.
    blogs = Post.objects.filter(status=Post.OPEN, type=Post.BLOG,
                                creation_date__gt=start).select_related('author')
    blogs = blogs[:limit]

    # Total post count
    total_post_count = Post.objects.filter(status=Post.OPEN).count()

    # Total user count
    total_user_count = User.objects.filter().count()

    hard_worker = User.objects.filter(post__status=Post.OPEN, post__lastedit_date__gt=start)\
        .annotate(total=Count("post")).order_by('-total').select_related("profile")
    hard_worker = hard_worker[:limit]


    params = dict(
        site = site,
        top_posts=top_posts,
        upd_posts=upd_posts,
        blogs=blogs,
        total_post_count=total_post_count,
        total_user_count=total_user_count,
        start=start.strftime("%b %d, %Y"),
        hard_worker=hard_worker,
        days=days,
    )

    text_body = html_body = ''

    if text_tmpl:
        text_body = html.render(text_tmpl, **params)
    if html_tmpl:
        html_body = html.render(html_tmpl, **params)

    if verbosity>0:
        print text_body
        print html_body

    # params = dict(subject=subject, from_email=from_email, recipient_list=recipient_list)
    #send_mail(subject=subject, message=message, from_email=from_email, recipient_list=recipient_list)


class Command(BaseCommand):
    help = 'send digest emails'

    option_list = BaseCommand.option_list + (

        make_option('--test', dest='test', action='store_true', default=False,
                    help='sends a test email to the site admin with the digest email (default=%default)'),

        make_option('--send', dest='send', action='store_true', default=False,
                    help='sends the email to all subscribers. Both text and html templates '
                         'are sent in a single multi-part email. (default=%default)'),

        make_option('--days', dest='days', default=1, type=int, metavar='NUMBER',
                    help='how many days back to include into the digest (default=%default)'),

        make_option('--limit', dest='limit', default=10, type=int, metavar='NUMBER',
                    help='limit the number of posts to this value (default=%default)'),

        make_option('--text', dest='text_template', metavar='FILE', default="digest/daily_digest.txt",
                    help='text template to use (default=%default)'),

        make_option('--html', dest='html_template',metavar='FILE', default=None,
                    help='html template to use (default=%default)'),
    )

    def handle(self, *args, **options):

        days, send, test = options['days'], options['send'], options['test']
        verbosity = options['verbosity']
        text_tmpl, html_tmpl =  options['text_template'], options['html_template']

        daily_digest(days=days, send=send, text_tmpl=text_tmpl, html_tmpl=html_tmpl,
            verbosity=int(verbosity),)




