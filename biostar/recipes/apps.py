import os
import csv
from django.db.models.signals import post_migrate
from django.core.files import File
from django.apps import AppConfig
from django.conf import settings
from biostar.accounts.apps import init_app


class EngineConfig(AppConfig):
    name = 'biostar.recipes'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_app, sender=self)
        post_migrate.connect(init_snippets, sender=self)
        pass


def init_snippets(sender, **kwargs):
    from biostar.recipes.models import SnippetType, Snippet
    from biostar.accounts.models import User

    data_root = os.path.join(settings.BASE_DIR, 'initial')

    # Get csv file with initial code snippets.
    default_snippets = os.path.join(data_root, 'initial-snippets.csv')

    # Get the owner of these snippets to be a superuser.
    owner = User.objects.filter(is_superuser=True).first()

    stream = open(default_snippets, 'rU')

    for row in csv.DictReader(stream):
        # Parse row content
        uid = row.get("category_uid", '').strip()
        snippet_uid = row.get("snippet_uid", '').strip()
        name = row.get("category_name", '').strip()
        image = row.get("image", '').strip()
        label = row.get('label', '').strip()
        command = row.get('command', '').strip()

        # Check if this default snippet type already exists.
        cmd_type = SnippetType.objects.filter(uid=uid).first()
        if cmd_type is None:
            # Create the default snippet type
            cmd_type = SnippetType.objects.create(name=name, uid=uid, owner=owner, default=True)
            if image:
                stream = File(open(os.path.join(data_root, image), 'rb'))
                cmd_type.image = stream
                cmd_type.save()

        # Check if the actual code snippet already exists
        snippet = Snippet.objects.filter(uid=snippet_uid).first()
        if snippet:
            continue

        # Create code snippet.
        Snippet.objects.create(help_text=label, uid=snippet_uid, default=True,
                               command=command, type=cmd_type, owner=owner)
