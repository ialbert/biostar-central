from django.apps import AppConfig
from django.db.models.signals import post_migrate
from biostar3.forum import models
from django.contrib.auth.models import Permission
from django.conf import settings
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

def post_migrate_tasks(sender, **kwargs):
    # The permissions should already exist. See comment on create_permissions.

    can_moderate_user = Permission.objects.get(codename=models.MODERATE_USER_PERMISSION)
    can_moderate_post = Permission.objects.get(codename=models.MODERATE_POST_PERMISSION)
    can_ban_user = Permission.objects.get(codename=models.BAN_USER_PERMISSION)

    # Ensure admin and moderator groups exist and have the right permissions.
    admin_group, admin_created = models.Group.objects.get_or_create(name=models.ADMIN_GROUP_NAME)

    if admin_created:
        admin_group.permissions.add(can_moderate_user, can_moderate_post, can_ban_user)

    mod_group, mod_created = models.Group.objects.get_or_create(name=models.MODERATOR_GROUP_NAME)
    if mod_created:
        mod_group.permissions.add(can_moderate_user, can_moderate_post)


    # Create the default admin user.
    for name, email in settings.ADMINS:
        admin, created = models.User.objects.get_or_create(email=email)
        if created:
            admin.name = name
            admin.is_staff = admin.is_admin = admin.is_superuser = True
            admin.type = models.User.ADMIN
            admin.set_password(settings.SECRET_KEY)
            admin.groups.add(admin_group, mod_group)
            admin.save()
            logger.info("added admin user with email=%s, password=SECRET_KEY, name=%s" % (admin.email, admin.name))

    # Initialize