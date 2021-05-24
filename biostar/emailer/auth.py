import logging
from .models import EmailSubscription, EmailGroup

logger = logging.getLogger("engine")


def add_subscription(email, group, name=''):

    # Fetch the subscriptions if these may exists.
    query = EmailSubscription.objects.filter(group=group, email=email)

    # Drop subscription if it exists.
    query.delete()

    # Create the new subscription.
    EmailSubscription.objects.create(group=group, email=email)
