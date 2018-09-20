
import logging

from django.core.management.base import BaseCommand

from biostar.engine.models import Access
from biostar.message.models import Message
from biostar.forum.models import Post
from biostar.accounts.models import User

logger = logging.getLogger("engine")


class Command(BaseCommand):
    help = "Add users"

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        # Check users who have access to a project
        users = Access.objects.filter(access__gte=Access.NO_ACCESS).values_list("user", flat=True)

        # Users to deleted
        users = User.objects.exclude(pk__in=users)

        # Delete posts first
        Post.objects.all().delete()

        # Delete messages next
        Message.objects.all().delete()

        # Delete the users
        users.delete()


