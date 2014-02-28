from django.db import models
from django.conf import settings

# Create your models here.

class Badge(models.Model):
    BRONZE, SILVER, GOLD = range(3)
    CHOICES = ((BRONZE, 'Bronze'), (SILVER, 'Silver'), (GOLD, 'Gold'))

    name = models.CharField(max_length=50)
    description = models.CharField(max_length=200)
    type = models.IntegerField(choices=CHOICES, default=BRONZE)

    # Unique badges may be earned only once
    unique = models.BooleanField(default=False)

    # Secret badges are not listed on the badge list
    secret = models.BooleanField(default=False)

    # Total number of times awarded
    count  = models.IntegerField(default=0)

    def get_absolute_url(self):
        return "/badge/show/%s/" % self.id

    def __unicode__(self):
        return self.name

class Award(models.Model):
    '''
    A badge being awarded to a user.Cannot be ManyToManyField
    because some may be earned multiple times
    '''
    badge = models.ForeignKey(Badge)
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    date = models.DateTimeField(auto_now_add=True)