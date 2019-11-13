import os
import logging

import hjson
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.recipes.models import Project, Access, Analysis
from biostar.recipes import util, auth

logger = logging.getLogger("engine")

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


def strip_json(json_text):
    """
    Strip settings parameter in json_text to only contain execute options
    Deletes the 'settings' parameter if there are no execute options.
    """
    try:
        local_json = hjson.loads(json_text)
    except Exception as exep:
        logger.error(f'Error loading json text: {exep}')
        return

    # Fetch the execute options
    execute_options = local_json.get('settings', {}).get('execute', {})

    # Check to see if it is present
    if execute_options:
        # Strip run settings of every thing but execute options
        local_json['settings'] = dict(execute=execute_options)
    else:
        # NOTE: Delete 'settings' from json text
        local_json['settings'] = ''
        del local_json['settings']

    new_json = hjson.dumps(local_json)
    return new_json


@receiver(post_save, sender=Project)
def finalize_project(sender, instance, created, raw, update_fields, **kwargs):
    # Ensure a project has at least one recipe on creation.
    if created and not instance.analysis_set.exists():
        # Add starter hello world recipe to project.
        try:
            json_text = open(join(DATA_DIR, 'starter.hjson'), 'r').read()
            template = open(join(DATA_DIR, 'starter.sh'), 'r').read()
            image = os.path.join(DATA_DIR, 'starter.png')
            image_stream = open(image, 'rb')
        except Exception as exc:
            logger.error(f'{exc}')
            json_text = '{}'
            template = "echo 'Hello World'"
            image_stream = None

        name = 'First recipe'
        text = "This recipe was created automatically."

        # Create starter recipe.
        auth.create_analysis(project=instance, json_text=json_text, template=template,
                             name=name, text=text, stream=image_stream)


@receiver(post_save, sender=Analysis)
def finalize_recipe(sender, instance, created, raw, update_fields, **kwargs):

    # Strip json of 'settings' parameter
    instance.json_text = strip_json(instance.json_text)
    root_is_writable = auth.writeable_recipe(user=instance.lastedit_user, source=instance, project=instance.project)

    if instance.is_cloned:
        root = instance.root
        # Final check to see the clone's last edit user
        # has write access to the root
        if root_is_writable:
            # Update root with instance data.
            root.merge(instance, save=True)
        return

    if instance.is_root:
        # Update information of all children belonging to this root.
        instance.update_children()
