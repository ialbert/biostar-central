import uuid

from django.db import models

MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_TEMPLATE_LEN = 20 * MAX_TEXT_LEN


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


class EmailGroup(models.Model):
    """
    Represents an group of email addresses.
    """
    name = models.CharField(max_length=MAX_NAME_LEN)
    uid = models.CharField(max_length=32, unique=True)
    text = models.CharField(max_length=MAX_TEXT_LEN)
    html = models.CharField(max_length=MAX_TEXT_LEN)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        self.uid = self.uid or get_uuid(16)
        super(EmailGroup, self).save()


class EmailSubscription(models.Model):
    """
    Connects email groups to email addresses.
    """
    ACTIVE, DELETED, UNSUBSCRIBE = 1, 2, 3
    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (UNSUBSCRIBE, "Unsubscribed")]

    uid = models.CharField(max_length=32, unique=True)
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)
    email = models.CharField(max_length=MAX_NAME_LEN, default='', blank=False)
    group = models.ForeignKey(EmailGroup, on_delete=models.CASCADE)

    def __str__(self):
        return f"{self.email} | {self.group.name}"

    def save(self, *args, **kwargs):
        self.uid = self.uid or get_uuid(16)
        super(EmailSubscription, self).save()

    def active(self):
        return self.state == self.ACTIVE
