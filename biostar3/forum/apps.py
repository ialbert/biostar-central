from django.apps import AppConfig
from django.db.models.signals import post_migrate
from biostar3.forum import models
from biostar3.forum.models import User, Group, Post
from django.contrib.auth.models import Permission
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from django.contrib.sites.models import Site
from django.db import transaction

import logging

logger = logging.getLogger('biostar')

# This is essential to be here as it registers
# the permission creation signal before our post_migrate
# connector below.
from django.contrib.auth.management import create_permissions

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


def check_permission(user, perm_name, label=BiostarAppConfig.app_label):
    """
    Generic function to apply a permission on the current app.
    """
    perm_label = "%s.%s" % (label, perm_name)
    return user.has_perm(perm_label)


# Functions to check permissions of the user on an action.
authorize_post_mod = lambda user: check_permission(user, models.MODERATE_POST_PERMISSION)
authorize_user_mod = lambda user: check_permission(user, models.MODERATE_USER_PERMISSION)
authorize_user_ban = lambda user: check_permission(user, models.BAN_USER_PERMISSION)

def post_migrate_tasks(sender, **kwargs):
    """
    Sets up data post migration. The site will rely on data set up via this function.
    """

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

    # The permissions should already exist. See comment on create_permissions.
    can_moderate_user = Permission.objects.get(codename=models.MODERATE_USER_PERMISSION)
    can_moderate_post = Permission.objects.get(codename=models.MODERATE_POST_PERMISSION)
    can_ban_user = Permission.objects.get(codename=models.BAN_USER_PERMISSION)

    # Get the first admin user
    admin = models.User.objects.filter(type=User.ADMIN).first()
    if not admin:
        raise ImproperlyConfigured("settings must include ADMINS attribute.")

    # Create the default group.
    default_group, default_created = models.get_or_create_group(name=settings.DEFAULT_GROUP_NAME, user=admin)

    # Update all toplevel posts with no groups to have the default group.
    Post.objects.filter(type__in=Post.TOP_LEVEL, group=None).update(group=default_group)

    # Admin group and premissions
    admin_group, admin_created = models.get_or_create_group(name=models.ADMIN_GROUP_NAME, user=admin)
    admin_group.permissions.add(can_moderate_user, can_moderate_post, can_ban_user)

    # Moderator group and permissions
    mod_group, mod_created = models.get_or_create_group(name=models.MODERATOR_GROUP_NAME, user=admin)
    mod_group.permissions.add(can_moderate_user, can_moderate_post)

    # All admin users need to have admin group level permissions.
    for user in models.User.objects.filter(type=User.ADMIN):
        user.groups.add(admin_group, mod_group)

    # All moderator users need to have moderator level permissions.
    for user in models.User.objects.filter(type=User.MODERATOR):
        user.groups.add(mod_group)

    # Sets up the default domain
    site = Site.objects.get_current()
    if site.domain != settings.SITE_DOMAIN:
        site.name = settings.SITE_NAME
        site.domain = settings.SITE_DOMAIN
        site.save()
        logger.info("adding site=%s, name=%s, domain=%s" % (site.id, site.name, site.domain))

    # This is only needed when migrating the database.
    # Migrate tags if these exist.
    for post in Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val=''):
        tags = post.tag_val.split(",")
        post.tags.set(*tags)

    # Reset tag_val field. This attribute will be dropped on a second migration.
    Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(tag_val='').update(tag_val='')