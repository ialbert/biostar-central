import os
import hjson
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.recipes.models import Project, Access, Analysis
from biostar.recipes import util, auth


__CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

DATA_DIR = os.path.join(__CURRENT_DIR, 'recipes')


def join(*args):
    return os.path.join(*args)


@receiver(post_save, sender=Project)
def update_access(sender, instance, created, raw, update_fields, **kwargs):
    # Give the owner WRITE ACCESS if they do not have it.
    entry = Access.objects.filter(user=instance.owner, project=instance, access=Access.WRITE_ACCESS)
    if entry.first() is None:
        entry = Access.objects.create(user=instance.owner, project=instance, access=Access.WRITE_ACCESS)


@receiver(post_save, sender=Project)
def finalize_project(sender, instance, created, raw, update_fields, **kwargs):

    # Ensure a project has at least one recipe on creation.
    if created and not instance.analysis_set.exists():
        # Add starter hello world recipe to project.

        json_text = open(join(DATA_DIR, 'starter.hjson'), 'r').read()
        template = open(join(DATA_DIR, 'starter.sh'), 'r').read()

        name = 'Starter Recipe'
        text = "This starter recipe prints **Hello World!** to the console."

        image = os.path.join(DATA_DIR, 'starter.png')
        image_stream = open(image, 'rb')
        # Create starter recipe.
        auth.create_analysis(project=instance, json_text=json_text, template=template,
                             name=name, text=text, stream=image_stream)