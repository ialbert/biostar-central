
import logging

from django.core.management.base import BaseCommand
from biostar.accounts.models import Profile
from django.contrib.auth import get_user_model

User = get_user_model()
logger = logging.getLogger("engine")


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        # Get users from default database

        all_users1 = User.objects.using("default").all()
        all_users2 = User.objects.using("users").all()

        print(all_users1.count())
        print(all_users2.count())

        for user1 in all_users1:
            # Create user in second user
            user2 = User.objects.using("users").filter(email=user1.email).first()

            if user2:
                continue

            #TODO: check if password is correctly migrated.
            user2 = User(username=user1.username, email=user1.email, password=user1.password)
            user2.save(using="users")
            Profile.objects.using("users").filter(user=user2).update(uid=user1.profile.uid, name=user1.profile.name,
                              max_upload_size=user1.profile.max_upload_size, role=user1.profile.role,
                              last_login=user1.profile.last_login, new_messages=user1.profile.new_messages,
                              date_joined=user1.profile.date_joined, location=user1.profile.location,
                              website=user1.profile.website, scholar=user1.profile.scholar,
                              score=user1.profile.score, twitter=user1.profile.twitter,
                              my_tags=user1.profile.my_tags, text=user1.profile.text,
                              html=user1.profile.html, email_verified=user1.profile.email_verified,
                              notify=user1.profile.notify, message_prefs=user1.profile.message_prefs,
                              digest_prefs=user1.profile.digest_prefs, opt_in=user1.profile.opt_in)
            print(user1.pk)

        all_users1 = User.objects.using("default").all()
        all_users2 = User.objects.using("users").all()

        print(all_users1.count())
        print(all_users2.count())
        return