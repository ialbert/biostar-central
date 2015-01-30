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


def post_migrate_tasks(sender, **kwargs):
    # The permissions should already exist. See comment on create_permissions.
    can_moderate = Permission.objects.get(codename=models.MODERATE_POST_PERMISSION)

    # Ensure admin and moderator groups exist and have the right permissions.
    admin_group, admin_created = models.Group.objects.get_or_create(name=models.ADMIN_GROUP_NAME)
    if admin_created:
        admin_group.permissions.add(can_moderate)

    mod_group, mod_created = models.Group.objects.get_or_create(name=models.MODERATOR_GROUP_NAME)
    if mod_created:
        mod_group.permissions.add(can_moderate)

    # Create the default admin user.
    for name, email in settings.ADMINS:
        admin, created = models.User.objects.get_or_create(email=email)
        if created:
            admin.name = name
            admin.is_staff = admin.is_admin = True
            admin.type = models.User.ADMIN
            admin.set_password(settings.SECRET_KEY)
            admin.save()
            logger.info("added admin user with email=%s, password=SECRET_KEY, name=%s" % (admin.email, admin.name))

    # Initialize