import logging
from django.conf import settings
logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    @timer(30000)
    def send_digest_messages():
        "Send posts as a digest"
        from biostar.forum.models import Subscription
        from biostar.message.models import Message

        digest_subs = Subscription.objects.filter(type=Subscription.DIGEST_MESSAGES)

        return


    @timer(300)
    def send_local_messages():
        "Send post subs as a local message"
        from biostar.message.models import Message
        from biostar.forum.models import Subscription

        local_subs = Subscription.objects.filter(type=Subscription.LOCAL_MESSAGE)

        return


    @timer(30)
    def send_email_messages():
        "Send posts as a email messages"
        from biostar.message.models import Message
        from biostar.forum.models import Subscription

        email_subs = Subscription.objects.filter(type=Subscription.LOCAL_MESSAGE)

        return


except ModuleNotFoundError as exc:
    pass
