from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging, glob, json, html2text
from biostar3.forum.models import *
from django.db.models import F, Q


def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


def perform_migration():
    """
    Migrates data from 2 to 3.
    Must be run when migrating a database from 2 to 3.
    """

    # Get the first admin user
    admin = User.objects.get(email=settings.ADMINS[0][1])

    # Update user scores
    logger.info("update user score")
    User.objects.all().update(score=F('score') * 10)


    # Save all users to trigger their html method.
    logger.info("resaving all user profiles")
    for prof in Profile.objects.all().exclude(info__isnull=True):
        try:
            prof.info = html2text.html2text(prof.info)
            prof.save()
        except Exception as exc:
            logger.error("error parsing profile %s" % prof.id)

    # Migrate old style tags to new style tags.
    logger.info('migrating tags')
    for post in Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val=''):
        tags = post.tag_val.split(",")
        post.tags.set(*tags)
        PostSub.objects.create(post=post, user=post.author)

    # Reset tag_val field.
    logger.info('resetting tag_val')
    Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val='').update(tag_val='')

    # All existing post will have the uuid be equal to their primary keys.
    logger.info("Add 'uuid' field to posts")
    Post.objects.all().update(uuid=F('pk'))

class Command(BaseCommand):
    help = '''Migrates data from version 2 to 3'''

    option_list = BaseCommand.option_list + (
        make_option('--all', action='store_true', default=False,
                    help='migrates data from version 2 to 3'),
    )

    def handle(self, *args, **options):
        if options["all"]:
            perform_migration()

