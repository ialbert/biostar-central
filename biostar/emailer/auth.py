import logging
from .models import EmailGroup, EmailAddress

logger = logging.getLogger("engine")




def add_sub(email, group, name=None):

    name = name or email
    # do a subscribtion_set.create(group=group)

    print(group)
    1/0
    logger.info(f"Subscribed {email} to ({group}) mailing-list.")
    pass



