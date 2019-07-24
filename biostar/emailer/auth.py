import logging
from .models import EmailAddress, Subscription

logger = logging.getLogger("engine")


def add_subscription(email, group, name=''):

    # Get the address from the database.
    address = EmailAddress.objects.filter(email=email).first()
    if not address:
        address = EmailAddress.objects.create(name=name, email=email)

    # Fetch the subscriptions if these may exists.
    query = Subscription.objects.filter(group=group, address=address)

    # Drop subscription if it exists.
    query.delete()
    # Create the new subscription.
    Subscription.objects.create(group=group, address=address)
