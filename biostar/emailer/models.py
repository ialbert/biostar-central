
from django.db import models
from biostar.accounts.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.engine.util import get_uuid


MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_TEMPLATE_LEN = 20 * MAX_TEXT_LEN



class EmailGroup(models.Model):

    name = models.CharField(max_length=MAX_NAME_LEN)
    uid =  models.CharField(max_length=16, blank=True, unique=True)
    text = models.CharField(max_length=MAX_TEXT_LEN)
    html = models.CharField(max_length=MAX_TEXT_LEN)

    # Mailing template specific to every group
    #template = models.CharField(max_length=MAX_TEMPLATE_LEN)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        self.uid = self.uid or get_uuid(16)

        super(EmailGroup, self).save()


class EmailAddress(models.Model):

    ACTIVE, DELETED, INACTIVE, UNSUBSCRIBE = 1,2,3,4
    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (UNSUBSCRIBE, "Unsubscirbed")]

    # required email
    email = models.CharField(max_length=MAX_NAME_LEN, unique=True, blank=False)
    name = models.CharField(max_length=MAX_NAME_LEN)

    uid = models.CharField(max_length=16, blank=True, unique=True)
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        self.uid = self.uid or get_uuid(16)

        super(EmailAddress, self).save()


class Subscription(models.Model):

    ACTIVE, DELETED, INACTIVE, UNSUBSCRIBE = 1,2,3,4
    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (UNSUBSCRIBE, "Unsubscirbed")]

    uid = models.CharField(max_length=16, blank=True, unique=True)
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)
    address = models.ForeignKey(EmailAddress)
    group = models.ForeignKey(EmailGroup)


    def save(self, *args, **kwargs):
        self.uid = self.uid or get_uuid(16)

        super(Subscription, self).save()
