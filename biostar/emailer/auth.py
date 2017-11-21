import logging
from .models import EmailGroup, EmailAddress, Subscription

logger = logging.getLogger("engine")



def add_sub(email, group, name=None):

    name = name or email
    address = EmailAddress.objects.filter(email=email).first()
    if not address:
        address = EmailAddress(name=name, email=email)
        address.save()

    sub = address.subscription_set.filter(group=group).first()

    if not sub:
        sub = address.subscription_set.create(group=group)
        logger.info(f"Subscribed {email} to ({group}) mailing-list.")
    else:
        logger.info(f"{email} already subscribed to ({group}) mailing-list.")

    return sub



