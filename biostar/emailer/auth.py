import logging
from .models import EmailGroup, EmailAddress

logger = logging.getLogger("engine")



def add_sub(email, group, name=None):

    name = name or email
    address = EmailAddress(name=name, email=email)
    address.save()
    sub = address.subscription_set.create(group=group)

    logger.info(f"Subscribed {email} to ({group}) mailing-list.")
    return sub



