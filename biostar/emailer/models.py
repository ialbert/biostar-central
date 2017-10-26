
from django.db import models
from engine.util import get_uuid



class EmailList(models.Model):

    ACTIVE, DELETED, NOT_ACTIVE = 1,2,3

    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (NOT_ACTIVE, "Not active")]

    uid = models.CharField(max_length=16, blank=True, unique=True, default=get_uuid(16))


