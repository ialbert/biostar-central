from django.apps import AppConfig
from django.db.models.signals import post_migrate
from biostar3.forum import models
from biostar3.forum.models import User, UserGroup, GroupPerm, Post
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from django.contrib.sites.models import Site
from django.db import transaction
from django.contrib.staticfiles import finders

import logging

logger = logging.getLogger('biostar')

class BiostarAppConfig(AppConfig):
    name = 'biostar3.forum'
    verbose_name = 'Biostar Forum'

    def ready(self):
        """
        Triggers when the application configuration is ready.
        """
        post_migrate.connect(post_migrate_tasks, sender=self)

    @property
    def app_label(self):
        return self.name.split(".")[-1]

def post_migrate_tasks(sender, **kwargs):
    """
    Sets up data post migration. The site will rely on data set up via this function.
    """

    # Check that the default logo exists.
    default_logo = finders.find(settings.DEFAULT_GROUP_LOGO)
    if not default_logo:
        raise ImproperlyConfigured("Cannot find default group logo at %s" % settings.DEFAULT_GROUP_LOGO)


    # Chicken and egg problem. Default group needs to be created before the first user.
    default_group, default_flag = UserGroup.objects.get_or_create(
        name=settings.DEFAULT_GROUP_NAME)



    # Create the default admin user.
    for name, email in settings.ADMINS:
        admin, created = User.objects.get_or_create(email=email)
        if created:
            admin.name = name
            admin.is_staff = admin.is_admin = admin.is_superuser = True
            admin.type = User.ADMIN
            admin.set_password(settings.SECRET_KEY)
            admin.save()
            logger.info("added admin user with email=%s, password=SECRET_KEY, name=%s" % (admin.email, admin.name))

    # Get the first admin user
    admin = models.User.objects.filter(type=User.ADMIN).first()
    if not admin:
        raise ImproperlyConfigured("settings must include ADMINS attribute.")

    # Associate the default group with the site admin.
    if default_flag:
        default_group.owner = admin
        default_group.save()

    # Create the meta group, talk about the site
    meta_name, meta_domain, meta_description = "Meta Group", "meta", "Discussions about the site itself"
    meta_group, meta_flag = UserGroup.objects.get_or_create(domain=meta_domain)
    if meta_flag:
        meta_group.name = meta_name
        meta_group.description = meta_description
        meta_group.owner = admin
        meta_group.save()

    # Update all toplevel posts with no groups to have the default group.
    logger.info("adding groups to posts")
    Post.objects.filter(type__in=Post.TOP_LEVEL, group=None).update(group=default_group)

    # All admin users need to have admin group level permissions.
    for user in models.User.objects.filter(type=User.ADMIN):
        GroupPerm.objects.get_or_create(group=default_group, user=user, type=GroupPerm.ADMIN)

    # All moderator users need to have moderator level permissions.
    for user in models.User.objects.filter(type=User.MODERATOR):
        GroupPerm.objects.get_or_create(group=default_group, user=user, type=GroupPerm.MODERATE)

    logger.info("adding groups to users")
    # Add group info to every user.
    for user in models.User.objects.all():
        user.usergroups.add(default_group)

    # Sets up the default domain
    site = Site.objects.get_current()
    if site.domain != settings.SITE_DOMAIN:
        site.name = settings.SITE_NAME
        site.domain = settings.SITE_DOMAIN
        site.save()
        logger.info("adding site=%s, name=%s, domain=%s" % (site.id, site.name, site.domain))

    # This is only needed when migrating the database.
    # Migrate tags if these exist.
    logger.info('migrating tags')
    for post in Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val=''):
        tags = post.tag_val.split(",")
        post.tags.set(*tags)

    # Reset tag_val field. This attribute will be dropped on a second migration.
    logger.info('resetting tag_val')
    Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val='').update(tag_val='')