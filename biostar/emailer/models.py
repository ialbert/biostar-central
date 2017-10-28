
from django.db import models
from biostar.engine.util import get_uuid



class EmailList(models.Model):

    ACTIVE, DELETED, INACTIVE, OPTOUT = 1,2,3,4

    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (OPTOUT, "Opted out")]

    uid = models.CharField(max_length=16, blank=True, unique=True, default=get_uuid(16))
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)

