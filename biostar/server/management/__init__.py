from django.db.models.signals import post_syncdb
from .commands.initialize_site import init_admin, init_domain, init_flatpages, init_social_providers

from south.signals import post_migrate

from biostar.server import models

# Populate the database with default admin, site info and social providers.
# This needs to run as early as possible to correctly populate the database.
def initialize(sender, **kwargs):
    init_admin()
    init_domain()
    init_social_providers()
    init_flatpages()

#post_syncdb.connect(initialize, sender=models)
