from django.db.models.signals import post_syncdb, pre_save
from biostar.apps.accounts.management.signals import initialize
from biostar.apps.accounts import models

# Populate the database with default admin, site info and social providers.
post_syncdb.connect(initialize, sender=models)
