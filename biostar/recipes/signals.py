import os
import logging

import mistune
import toml
from django.conf import settings
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


def strip_json(json_text):
    """
    Strip settings parameter in json_text to only contain execute options
    Deletes the 'settings' parameter if there are no execute options.
    """

    try:
        local_dict = toml.loads(json_text)
    except Exception as exep:
        logger.error(f'Error loading json text: {exep}, {json_text}')
        return json_text

    # Fetch the execute options
    execute_options = local_dict.get('settings', {}).get('execute', {})
    data_options = local_dict.get('settings', {}).get('create', {})

    # Check to see if it is present
    if execute_options or data_options:
        # Strip run settings of every thing but execute options
        local_dict['settings'] = dict(execute=execute_options, create=data_options)
    else:
        # NOTE: Delete 'settings' from json text
        local_dict['settings'] = ''
        del local_dict['settings']

    new_json = toml.dumps(local_dict)

    return new_json


def initial_recipe(project):
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
    recipe = auth.create_analysis(project=project, json_text=json_text, template=template,
                                  name=name, text=text, stream=image_stream, security=Analysis.AUTHORIZED)
    return recipe


@receiver(post_save, sender=Project)
def finalize_project(sender, instance, created, raw, update_fields, **kwargs):

    if created:
        # Generate friendly uid
        uid = auth.new_uid(obj=instance, objtype=Project, prefix="project")
        instance.uid = uid

        # Set the project directory
        instance.dir = instance.dir or join(settings.MEDIA_ROOT, "projects", f"{instance.uid}")

        # Get project with highest rank and add to it,
        # ensuring this new project is at the top of lists
        first = Project.objects.order_by('-rank').first()
        instance.rank = first.rank + instance.pk if first else instance.pk

        # Create the job directory if it does not exist.
        os.makedirs(instance.dir, exist_ok=True)

        # Update project fields.
        Project.objects.filter(id=instance.id).update(uid=instance.uid, dir=instance.dir, rank=instance.rank)

        # Create a starter recipe if none exist.
        if not instance.analysis_set.exists():
            initial_recipe(project=instance)

    # Cascade deleted states to recipe, data, and results.
    if instance.deleted:
        Analysis.objects.filter(project__id=instance.pk).update(deleted=True)
        Data.objects.filter(project__id=instance.pk).update(deleted=True)
        Job.objects.filter(project__id=instance.pk).update(deleted=True)


@receiver(post_save, sender=Analysis)
def finalize_recipe(sender, instance, created, raw, update_fields, **kwargs):

    if created:
        # Generate friendly uid
        uid = auth.new_uid(obj=instance, objtype=Analysis, prefix="recipe")
        instance.uid = uid

        Analysis.objects.filter(id=instance.id).update(uid=instance.uid)

        # Get recipe with highest rank and add to it,
        # ensuring this new recipe is at the top of lists
        first = Analysis.objects.order_by('-rank').first()
        instance.rank = first.rank + instance.pk if first else instance.pk

    # Update the last edit date and user of project
    user = instance.lastedit_user

    # Strip json text of 'settings' parameter
    instance.json_text = strip_json(instance.json_text)
    Project.objects.filter(id=instance.project.id).update(lastedit_date=instance.lastedit_date,
                                                          lastedit_user=user)

    # Update information of all children belonging to this root.
    if instance.is_root:
        instance.update_children()

    # Update the project count and last edit date when job is created
    instance.project.set_counts()


@receiver(post_save, sender=Job)
def finalize_job(sender, instance, created, raw, update_fields, **kwargs):

    # Update the project count.
    instance.project.set_counts()

    if created:
        # Generate friendly uid
        uid = auth.new_uid(obj=instance, objtype=Job, prefix="job")
        instance.uid = uid
        # Generate the path based on the uid.
        instance.path = join(settings.MEDIA_ROOT, "jobs", f"{instance.uid}")

        # Create the job directory if it does not exist.
        os.makedirs(instance.path, exist_ok=True)

        # Update the information in db.
        Job.objects.filter(id=instance.id).update(uid=instance.uid, path=instance.path,
                                                  text=instance.text, html=instance.html)


@receiver(post_save, sender=Data)
def finalize_data(sender, instance, created, raw, update_fields, **kwargs):

    # Update the projects last edit user when a data is uploaded
    Project.objects.filter(id=instance.project.id).update(lastedit_user=instance.lastedit_user,
                                                          lastedit_date=instance.lastedit_date)
    # Update the project count.
    instance.project.set_counts()

    if created:
        # Generate friendly uid
        uid = auth.new_uid(obj=instance, objtype=Data, prefix="data")
        instance.uid = uid
        
        # Set the data directory with the recently created uid
        instance.dir = join(instance.get_project_dir(), f"{instance.uid}")

        # Set the toc file with the recently created uid
        instance.toc = join(settings.TOC_ROOT, f"toc-{instance.uid}.txt")

        # Build the data directory.
        os.makedirs(instance.dir, exist_ok=True)

        # Set the table of contents for the data
        if not os.path.isfile(instance.toc):
            with open(instance.toc, 'wt') as fp:
                pass

        # Update the dir, toc, and uid.
        Data.objects.filter(id=instance.id).update(uid=instance.uid, dir=instance.dir, toc=instance.toc)

    instance.make_toc()
