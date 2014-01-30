from django.db.models.signals import post_syncdb, pre_save
from biostar.apps.people.management import actions
from biostar.apps.people import models

# Populate the database with default admin, site info and social providers.
# This needs to run as early as possible to correctly populate the database.
post_syncdb.connect(actions.initialize, sender=models)
