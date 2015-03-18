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
from biostar3.forum import awards

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

            for award in awards.get_awards():
                logger.info("award %s" % award.uuid)
                award.check()