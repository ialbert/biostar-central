
import logging
from django.core.management.base import BaseCommand
from biostar.accounts.models import User, Group
from biostar.accounts import util

logger = logging.getLogger("engine")

class Command(BaseCommand):
    help = "Add users"
    def add_arguments(self, parser):

        parser.add_argument('--email', required=True,
                            help="Email of the user to be added")

        parser.add_argument('--username',
                            help="Username of the user to be added. defaults to random")

        parser.add_argument('--password', default="testbuddy11",
                            help="password associated with the email")

        parser.add_argument('--name',
                            help="User first name, defaults to email prefix.")

        parser.add_argument('--group_name',
                            help="User first name, defaults to email prefix.")

        parser.add_argument('--group_id',
                            help="User first name, defaults to email prefix.")

        parser.add_argument('--is_superuser',action='store_true', default=False,
                            help="Created user is an admin")

        parser.add_argument('--is_staff', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")


    def handle(self, *args, **options):

        email = options['email']
        username = options['username'] or util.get_uuid(8)
        password = options['password']
        name = options["name"] or email.split("@")[0]
        is_superuser = options["is_superuser"]
        is_staff = options["is_staff"]
        group_name = options["group_name"]
        group_id = options["group_id"]

        group = Group.objects.filter()
        if User.objects.filter(email=email).exists():
            logger.info(f"User with email={email} already exists")
            return
        user = User.objects.create(email=email, username=username, first_name=name,
                                   is_staff=is_staff, is_superuser=is_superuser)
        user.set_password(password)
        user.save()
        logger.info(f"Created user.email={email}, user.id={user.id}.")

