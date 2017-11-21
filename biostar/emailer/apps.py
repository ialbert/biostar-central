
from django.db.models.signals import post_migrate
from django.apps import AppConfig


class EmailerConfig(AppConfig):
    name = 'biostar.emailer'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_group, sender=self)
        post_migrate.connect(init_sub, sender=self)


def init_group(sender, **kwargs):
    "Initialize with 'staff' group/mailing-list"

    from biostar.emailer.models import EmailGroup

    if not EmailGroup.objects.filter(name="staff"):
        staff_group = EmailGroup(name="staff",text="Mailing list for staff")
        staff_group.save()
    return


def init_sub(sender, **kwargs):
    "Make a subscription to staff mailing list as a test"

    from biostar.emailer.models import EmailGroup, EmailAddress

    address = EmailAddress.objects.filter(email="testbuddy@lvh.me").first()
    if not address:
        address = EmailAddress(email="testbuddy@lvh.me", name="testbud")
        address.save()

        # Create a new subscription via the now created address object
        mailing_list = EmailGroup.objects.filter(name="staff").first()
        sub = address.subscription_set.create(group=mailing_list)

    return