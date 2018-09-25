import logging

from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger('engine')


class MessageConfig(AppConfig):
    name = 'biostar.message'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_messages, sender=self)

        pass


def init_messages(sender, **kwargs):

    from django.contrib.auth import get_user_model
    from . import auth, models

    User = get_user_model()

    name, email = settings.ADMINS[0]

    sender = User.objects.filter(email=email).first()
    body = "Hello from the biostar-engine developers, we hope you enjoy the website."
    subject = "Welcome to the biostar-engine!"

    test_2 = User.objects.filter(username="testing").first()
    if not test_2:
        # Create user and send message once.
        test_2 = User.objects.create(username="testing", email="testing@test")
        recipient_list = [sender, test_2]
        msg = auth.create_messages(body=body, subject=subject, recipient_list=recipient_list,
                            sender=sender)

        # Test with a message tree whenever debugging
        if settings.DEBUG:
            msg1 = msg[1]
            msg2 = msg[0]
            models.Message.objects.filter(pk=msg2.pk).update(parent_msg=msg1)
