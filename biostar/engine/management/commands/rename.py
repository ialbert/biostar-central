import os
from biostar.engine.models import Data, Project, Analysis
from django.core.management.base import BaseCommand


def correct_image_path(single_obj=None, queryset=[]):
    """Make sure the object's image path matches the project directory."""

    if single_obj and single_obj.image:
        single_obj.image = str(single_obj.image).replace("proj-", "")
        single_obj.save()

    for obj in queryset:
        if obj.image:
            obj.image = str(obj.image).replace("proj-", "")
            obj.save()


def rename(entry, replacing=""):

    new_name = entry.name.replace(replacing, "")
    new_path = os.path.join(os.path.dirname(entry.path), new_name)
    print(entry.name, new_name, new_path)
    try:
        os.rename(entry.path, new_path)
    except Exception as exc:
        print(f"Renaming error: {exc}")

    return new_path


class Command(BaseCommand):
    help = 'Make corrections without migrating'

    def add_arguments(self, parser):
        parser.add_argument('--root', required=True, help="Root directory with folders that need to be renamed.")

    def handle(self, *args, **options):
        root = options['root']

        for entry in os.scandir(root):
            if entry.name.startswith("job-"):
                rename(entry=entry, replacing="job-")

            elif entry.name.startswith("proj-"):

                proj_path = rename(entry=entry, replacing="proj-")
                proj_uid = os.path.basename(proj_path)
                project = Project.objects.filter(uid=proj_uid).first()

                correct_image_path(single_obj=project)

                # Go one level deeper in project directories to rename the data folders.
                for data in os.scandir(proj_path):
                    data_path = rename(entry=data, replacing="store-")
                    data_uid = os.path.basename(data_path)
                    data = Data.objects.filter(uid=data_uid).first()
                    if data:
                        # Data objects need their table of contents remade.
                        data.make_toc()

        recipes = Analysis.objects.get_all()

        # Recipes need to have image field with correct path
        correct_image_path(queryset=recipes)