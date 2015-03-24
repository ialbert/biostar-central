from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging, glob, json, html2text
from biostar3.forum.models import *


def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


def perform_migration():
    """
    Migrates data from 2 to 3.
    Must be run when migrating a database from 2 to 3.
    """

    # Get the default group.
    default_group = UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

    # Get the first admin user
    admin = User.objects.get(email=settings.ADMINS[0][1])

    # Create the meta group, talk about the site
    logger.info("creating meta group for admin")
    meta_name, meta_domain, meta_description = "Meta Talk", "meta", "Discussions about the site itself"
    meta_group, meta_flag = UserGroup.objects.get_or_create(domain=meta_domain)

    if meta_flag:
        meta_group.name = meta_name
        meta_group.description = meta_description
        meta_group.owner = admin
        meta_group.save()

    # Update all toplevel posts with no groups to have the default group.
    logger.info("adding groups to posts")
    Post.objects.filter(type__in=Post.TOP_LEVEL, usergroup=None).update(usergroup=default_group)

    # All admin users belong to the meta group and need to have admin group level permissions.
    for user in User.objects.filter(type=User.ADMIN).exclude(pk=admin.id):
        GroupSub.objects.get_or_create(user=user, usergroup=meta_group)
        GroupPerm.objects.get_or_create(usergroup=default_group, user=user, role=GroupPerm.ADMIN)

        # All moderator users need to have moderator level permissions.
    for user in User.objects.filter(type=User.MODERATOR).exclude(pk=admin.id):
        GroupSub.objects.get_or_create(user=user, usergroup=meta_group)
        GroupPerm.objects.get_or_create(usergroup=default_group, user=user, role=GroupPerm.MODERATE)

    # Add group info to every user.
    logger.info("adding groups to users")

    def group_generator():
        for user in User.objects.all().exclude(pk=admin.id):
            yield GroupSub(user=user, usergroup=default_group)

    # Bulk insert the groups.
    GroupSub.objects.bulk_create(group_generator())

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

    # Add group to all blogs that don't have one.
    Blog.objects.filter(usergroup=None).update(usergroup=default_group)

class Command(BaseCommand):
    help = '''Migrates data from version 2 to 3'''

    option_list = BaseCommand.option_list + (
        make_option('--all', action='store_true', default=False,
                    help='migrates data from version 2 to 3'),
    )

    def handle(self, *args, **options):
        if options["all"]:
            perform_migration()

