from biostar.engine.models import Data, Project, Analysis
from django.core.management.base import BaseCommand


def correct_image_path(queryset):
    """Make sure the object's image path matches the project directory."""

    for obj in queryset:
        if obj.image:
            obj.image = str(obj.image).replace("proj-", "")
            obj.save()


class Command(BaseCommand):
    help = 'Make corrections without migrating'

    def add_arguments(self, parser):

        pass

    def handle(self, *args, **options):

        projects = Project.objects.get_all()
        recipes = Analysis.objects.get_all()
        data = Data.objects.get_all()

        # Recipes and projects need to have image field with correct path
        correct_image_path(queryset=recipes)
        correct_image_path(queryset=projects)

        # Data objects need their table of contents remade.
        for obj in data:
            obj.make_toc()



