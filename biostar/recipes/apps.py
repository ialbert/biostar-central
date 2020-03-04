import os
import csv
from django.db.models.signals import post_migrate
#from django.contrib.redirects.models import Redirect
from django.core.files import File
from django.apps import AppConfig
from django.conf import settings
from biostar.accounts.apps import init_app


def join(*args):
    return os.path.abspath(os.path.join(*args))



class EngineConfig(AppConfig):
    name = 'biostar.recipes'

    def ready(self):
        from . import signals
        # Triggered upon app initialization.
        post_migrate.connect(init_app, sender=self)
        post_migrate.connect(init_dirs, sender=self)


def init_dirs(sender, **kwargs):
    from biostar.recipes.models import Project, Data

    projects = Project.objects.all()
    for project in projects:
        project.directory = join(settings.MEDIA_ROOT, "projects", f"{project.uid}")
        project.save()

    data = Data.objects.all()
    for datum in data:
        datum.directory = join(datum.project.dir, f"{datum.uid}")
        datum.toc = join(settings.TOC_ROOT, f"toc-{datum.uid}.txt")
        data.save()

        #Data.objects.filter(id=datum.project.id).update(dir=directory, toc=toc)
