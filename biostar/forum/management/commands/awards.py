import logging
from datetime import timedelta
from django.utils import timezone
from django.core.management.base import BaseCommand

from biostar.accounts.models import User
from biostar.forum.models import Award, Badge, Post
from biostar.forum.awards import ALL_AWARDS


logger = logging.getLogger('engine')


def give_award(user):

    # Get all valid awards.

    for award in ALL_AWARDS:
        # Valid award targets the user has earned
        targets = award.validate(user)

        for target in targets:
            post = target if isinstance(target, Post) else None
            badge = Badge.objects.filter(name=award.name).first()

            # Do not award a post multiple times.
            already_awarded = Award.objects.filter(user=user, badge=badge, post=post).exists()
            if post and already_awarded:
                continue

            # Create an award for each target.
            Award.objects.create(user=user, badge=badge, post=post)

            logger.info(f"award {badge.name} created for {user.email}" )



def create_user_awards(clear=False):

    # Clear all awards.
    if clear:
        Award.objects.all().delete()

    # GGet users active within the last day.
    past_day = timezone.now() - timedelta(days=1)

    users = User.objects.filter(profile__last_login__gte=past_day)
    logger.info(f"Checking awards for {users.count()} users.")

    for user in users:
        give_award(user)

    return


class Command(BaseCommand):
    help = 'Give awards to users.'

    def add_arguments(self, parser):

        parser.add_argument('-a', '--award_users', action='store_true', default=False,
                            help="Hand out awards to users.")
        parser.add_argument('-c', '--clear', action='store_true', default=False,
                            help="Clear current awards")

    def handle(self, *args, **options):

        give_awards = options['award_users']
        clear = options['clear']

        if give_awards:
            create_user_awards(clear=clear)
