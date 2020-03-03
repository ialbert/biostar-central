import os
import logging

import toml
import hjson
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.recipes.models import Project, Access, Analysis, Job, Data
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


def load_text(text):
    """
    Load text into a data dict.
    """
    try:
        # Try and load text as toml
        data = toml.loads(text)
    except Exception:
        # Try to load text as json
        data = hjson.loads(text)

    return data


def strip_json(json_text):
    """
    Strip settings parameter in json_text to only contain execute options
    Deletes the 'settings' parameter if there are no execute options.
    """
    try:
        local_dict = load_text(json_text)
    except Exception as exep:
        logger.error(f'Error loading json text: {exep}')
        return ""

    # Fetch the execute options
    execute_options = local_dict.get('settings', {}).get('execute', {})

    # Check to see if it is present
    if execute_options:
        # Strip run settings of every thing but execute options
        local_dict['settings'] = dict(execute=execute_options)
    else:
        # NOTE: Delete 'settings' from json text
        local_dict['settings'] = ''
        del local_dict['settings']

    new_json = toml.dumps(local_dict)
    return new_json


@receiver(post_save, sender=Project)
def finalize_project(sender, instance, created, raw, update_fields, **kwargs):

    # Ensure a project has at least one recipe on creation.
    if created and not instance.analysis_set.exists():
        # Generate friendly uid
        uid = auth.generate_uuid(prefix="project", suffix=instance.id)
        instance.uid = uid
        instance.label = uid
        Project.objects.filter(id=instance.id).update(uid=instance.uid, label=instance.label)

        # Add starter hello world recipe to project.
        try:
            json_text = open(join(DATA_DIR, 'starter.hjson'), 'r').read()
            template = open(join(DATA_DIR, 'starter.sh'), 'r').read()
            image = os.path.join(DATA_DIR, 'starter.png')
            image_stream = open(image, 'rb')
        except Exception as exc:
            logger.error(f'{exc}')
            json_text = ''
            template = "echo 'Hello World'"
            image_stream = None

        name = 'First recipe'
        text = "This recipe was created automatically."

        # Create starter recipe.
        auth.create_analysis(project=instance, json_text=json_text, template=template,
                             name=name, text=text, stream=image_stream)


@receiver(post_save, sender=Analysis)
def finalize_recipe(sender, instance, created, raw, update_fields, **kwargs):
    # Generate friendly uid
    if created:
        instance.uid = auth.generate_uuid(prefix="recipe", suffix=instance.id)
        Analysis.objects.filter(id=instance.id).update(uid=instance.uid)

    # Strip json of 'settings' parameter
    instance.json_text = strip_json(instance.json_text)
    # Update information of all children belonging to this root.
    if instance.is_root:
        instance.update_children()


@receiver(post_save, sender=Job)
def finalize_job(sender, instance, created, raw, update_fields, **kwargs):
    # Generate friendly uid
    if created:
        instance.uid = auth.generate_uuid(prefix="job", suffix=instance.id)
        Job.objects.filter(id=instance.id).update(uid=instance.uid)

        # Update the count and last edit date when job is created
        job_count = Job.objects.filter(deleted=False, project=instance.project).count()
        Project.objects.filter(id=instance.project.id).update(lastedit_user=instance.owner,
                                                              lastedit_date=util.now(),
                                                              jobs_count=job_count)


@receiver(post_save, sender=Data)
def finalize_data(sender, instance, created, raw, update_fields, **kwargs):
    # Generate friendly uid
    if created:
        instance.uid = auth.generate_uuid(prefix="data", suffix=instance.id)
        Data.objects.filter(id=instance.id).update(uid=instance.uid)
