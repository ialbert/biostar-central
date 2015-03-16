from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging, glob, json
import html2text
from biostar3.forum.models import Post
from biostar3.forum.html import sanitize

"""
Saves html to markdown format.
"""

logger = logging.getLogger('biostar')

def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = '''Parse html content into markdown.'''

    option_list = BaseCommand.option_list + (
        make_option('--all', action='store_true', default=False,
                    help='deletes sqlite database'),
    )

    def handle(self, *args, **options):

        if options['all']:

            post_count = Post.objects.all().count()

            logger.info("Resaving %s posts" % post_count)
            last = current = -1
            for index, post in enumerate(Post.objects.all().select_related("author")):
                try:
                    content = html2text.html2text(post.content)
                    html = sanitize(content, user=post.author)
                    Post.objects.filter(pk=post.id).update(content=content)
                    current = int(100.0 * index/post_count)
                    if last != current:
                        last = current
                        print ("progress: %d%%" % current)

                except Exception, exc:
                    logger.error(exc)
                    logger.error("Error parsing post %s" % post.id)