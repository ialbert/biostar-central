import os
from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings


class PlanetConfig(AppConfig):
    name = 'biostar.planet'

    def ready(self):
        #from . import signals
        # Triggered upon app initialization.
        #post_migrate.connect(init_awards, sender=self)
        #post_migrate.connect(init_digest, sender=self)
        post_migrate.connect(init_planet, sender=self)


def init_planet(sender, **kwargs):
    from biostar.planet import auth, models

    #models.Blog.objects.all().delete()
    #models.BlogPost.objects.all().delete()

    if settings.INIT_PLANET:
        # Load an initial set of blog posts
        fname = os.path.abspath(os.path.join(settings.PLANET_DIR, 'example-feeds.txt'))

        # Add feed to database if they do not exist already
        auth.add_blogs(fname)
