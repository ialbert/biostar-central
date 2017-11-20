
from django.db import models
from biostar.accounts.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.engine.util import get_uuid

MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000



class EmailGroup(models.Model):

    name = models.CharField(max_length=MAX_NAME_LEN)
    uid =  models.CharField(max_length=16, blank=True, unique=True, default=get_uuid(16))
    text = models.CharField(max_length=MAX_TEXT_LEN)
    html = models.CharField(max_length=MAX_TEXT_LEN)


class EmailAddress(models.Model):

    ACTIVE, DELETED, INACTIVE, Unsubscirbed = 1,2,3,4

    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (Unsubscirbed, "Unsubscirbed")]

    # required email
    email = models.CharField(max_length=MAX_NAME_LEN, unique=True, blank=False)

    name = models.CharField(max_length=MAX_NAME_LEN)

    uid = models.CharField(max_length=16, blank=True, unique=True, default=get_uuid(16))
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)


class Subscription(models.Model):

    ACTIVE, DELETED, INACTIVE, Unsubscirbed = 1,2,3,4

    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (Unsubscirbed, "Unsubscirbed")]

    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)
    address = models.ForeignKey(EmailAddress)
    group = models.ForeignKey(EmailGroup)
