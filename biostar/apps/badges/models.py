from django.db import models
from django.conf import settings
from django.core.urlresolvers import reverse
from django.core import mail
import logging

logger = logging.getLogger(__name__)

# Create your models here.

class Badge(models.Model):
    BRONZE, SILVER, GOLD = range(3)
    CHOICES = ((BRONZE, 'Bronze'), (SILVER, 'Silver'), (GOLD, 'Gold'))

    # The name of the badge.
    name = models.CharField(max_length=50)

    # The description of the badge.
    desc = models.CharField(max_length=200, default='')

    # The rarity of the badge.
    type = models.IntegerField(choices=CHOICES, default=BRONZE)

    # Unique badges may be earned only once
    unique = models.BooleanField(default=False)

    # Total number of times awarded
    count = models.IntegerField(default=0)

    # The icon to display for the badge.
    icon = models.CharField(default='fa fa-asterisk', max_length=250)

    def get_absolute_url(self):
        url = reverse("badge-details", kwargs=dict(pk=self.id))
        return url

    def __unicode__(self):
        return self.name


class Award(models.Model):
    '''
    A badge being awarded to a user.Cannot be ManyToManyField
    because some may be earned multiple times
    '''
    badge = models.ForeignKey(Badge)
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    date = models.DateTimeField()
    context = models.CharField(max_length=1000, default='')

class AwardDef(object):
    def __init__(self, name, desc, func, icon, type=Badge.BRONZE):
        self.name = name
        self.desc = desc
        self.fun = func
        self.icon = icon
        self.template = "badge/default.html"
        self.type = type

    def validate(self, *args, **kwargs):
        try:
            value = self.fun(*args, **kwargs)
            return value
        except Exception, exc:
            logger.error("validator error %s" % exc)
        return 0

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)


