from __future__ import unicode_literals, absolute_import, division
import logging, string, time

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
from django.core.urlresolvers import reverse
from django.core.mail import EmailMultiAlternatives

def render_digest(days, text_tmpl, html_tmpl, send, options, limit=10, verbosity=1):
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

    hard_worker = User.objects.filter(post__status=Post.OPEN, post__lastedit_date__gt=start) \
        .annotate(total=Count("post")).order_by('-total').select_related("profile")
    hard_worker = hard_worker[:limit]

    params = dict(
        site=site,
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

    if verbosity > 0:
        extras = dict(
            digest_manage=reverse("digest_manage"),
            digest_unsubscribe=reverse("digest_unsubscribe", kwargs=dict(uuid=1))
        )
        print text_body % extras
        print html_body % extras

    if send:

        logger.info('sending emails')
        emails = map(string.strip, open(send))

        def chunks(data, size):
            "Break into chunks of 100"
            for i in xrange(0, len(data), size):
                yield data[i:i+size]

        for chunk in chunks(emails, 100):
            users = User.objects.filter(email__in=chunk).select_related('profile')
            for user in users:
                try:
                    extras = dict(
                        digest_manage=reverse("digest_manage"),
                        digest_unsubscribe=reverse("digest_unsubscribe", kwargs=dict(uuid=user.profile.uuid))
                    )
                    text_content = text_body % extras
                    html_content = html_body % extras
                    subject = options['subject']
                    
                    msg = EmailMultiAlternatives(subject, text_content, from_email, [user.email])
                    msg.attach_alternative(html_content, "text/html")
                    msg.send()
                    time.sleep(0.3) # Throttle on Amazon.
                except Exception, exc:
                    logger.error('error %s sending email to %s' % (exc, user.email))


class Command(BaseCommand):
    help = 'send digest emails'

    option_list = BaseCommand.option_list + (

        make_option('--emails', dest='emails', default=None,
                    help='returns the emails for users with a certain digest preference set (daily, weekly, monthly)'),

        make_option('--subject', dest='subject', default="Biostar Digest",
                    help='the subject line of the email'),

        make_option('--send', dest='send', metavar="FILE", default=None,
                    help='sends the rendered templates to the list of emails in the file. Both text and html templates '
                         'are sent in a single multi-part email.'),

        make_option('--days', dest='days', default=1, type=int, metavar='NUMBER',
                    help='how many days back to include into the digest (default=%default)'),

        make_option('--limit', dest='limit', default=10, type=int, metavar='NUMBER',
                    help='limit the number of posts to this value (default=%default)'),

        make_option('--text', dest='text_template', metavar='FILE', default="digest/daily_digest.txt",
                    help='text template to use (default=%default)'),

        make_option('--html', dest='html_template', metavar='FILE', default=None,
                    help='html template to use (default=%default)'),

        make_option('--show', dest='show',
                    action="store_true", default=False, help='show the results'),

    )

    def handle(self, *args, **options):

        days, send = options['days'], options['send']
        emails, show = options['emails'], options['show']

        verbosity = options['verbosity']
        text_tmpl, html_tmpl = options['text_template'], options['html_template']

        if emails:
            # Shows the emails of users registered for a digest
            pass

        if show:
            # Produces the daily digest.
            render_digest(days=days, send=send, text_tmpl=text_tmpl, html_tmpl=html_tmpl, options=options,
                     verbosity=int(verbosity), )




