from django.apps import AppConfig
from django.db.models.signals import post_migrate
from biostar3.forum import models
from biostar3.forum.models import User, Group, Post
from django.contrib.auth.models import Permission
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured

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
        post_migrate.connect(post_migrate_tasks, sender=self)

    @property
    def app_label(self):
        return self.name.split(".")[-1]


def check_permission(user, perm_name):
    perm_label = "%s.%s" % (BiostarAppConfig.app_label, perm_name)
    return user.has_perm(perm_label)


authorize_post_mod = lambda user: check_permission(user, models.MODERATE_POST_PERMISSION)
authorize_user_mod = lambda user: check_permission(user, models.MODERATE_USER_PERMISSION)
authorize_user_ban = lambda user: check_permission(user, models.BAN_USER_PERMISSION)


def create_group(name, user):
    group, flag = models.Group.objects.get_or_create(name=name)
    if flag:
        models.GroupInfo.objects.create(group=group, author=user)
    return group, flag


def post_migrate_tasks(sender, **kwargs):
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

    # General group created.
    general_group, general_created = create_group(name=settings.DEFAULT_GROUP_NAME, user=admin)

    # Update all toplevel posts with no groups to have the general group.
    Post.objects.filter(type__in=Post.TOP_LEVEL, group=None).update(group=general_group)

    # Ensure admin and moderator groups exist and have the right permissions.
    admin_group, admin_created = create_group(name=models.ADMIN_GROUP_NAME, user=admin)
    admin_group.permissions.add(can_moderate_user, can_moderate_post, can_ban_user)

    mod_group, mod_created = create_group(name=models.MODERATOR_GROUP_NAME, user=admin)
    mod_group.permissions.add(can_moderate_user, can_moderate_post)

    # Update all admin users to have permissions.
    for user in models.User.objects.filter(type=User.ADMIN):
        user.groups.add(admin_group, mod_group)