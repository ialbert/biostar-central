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

    name = models.CharField(max_length=50)
    description = models.CharField(max_length=200)
    type = models.IntegerField(choices=CHOICES, default=BRONZE)

    # Unique badges may be earned only once
    unique = models.BooleanField(default=False)

    # Secret badges are not listed on the badge list
    secret = models.BooleanField(default=False)

    # Total number of times awarded
    count = models.IntegerField(default=0)

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
    date = models.DateTimeField(auto_now_add=True)

class BadgeDef(object):
    def __init__(self, name, desc, func):
        self.name = name
        self.desc = desc
        self.fun = func
        self.template = "badge/default.html"

    def validate(self, *args, **kwargs):
        try:
            return self.fun(*args, **kwargs)
        except Exception, exc:
            logger.error("validator error %s" % exc)

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)

# Simple badges can be computed via a single function call
AUTOBIO = BadgeDef(
    name = "Autobiographer",
    desc = "having more than 80 characters in the information field of the user's profile",
    func = lambda user: (len(user.profile.info) > 80)
)

SIMPLE_USER_BADGES = [
    AUTOBIO,
]

def init_badges():
    "Initializes the badges"
    logger.info("initializing badges")
    for obj in SIMPLE_USER_BADGES:
        badge = Badge.objects.get_or_create(name=obj.name, description=obj.desc)

# Tries to award a badge to the user
def create_user_award(user):
    "The callback is a function that called at the end of awardin"

    # Debug only
    #Award.objects.all().delete()

    # The awards the user has won at this point
    awards = set(a.badge.name for a in Award.objects.filter(user=user).select_related('badge'))

    # Function to detect if the award already exists
    not_awarded = lambda name: name not in awards

    for obj in SIMPLE_USER_BADGES:

        if not_awarded(obj.name) and obj.validate(user):

            # Get the badge
            badge = Badge.objects.get(name=obj.name)
            badge.count += 1
            badge.save()

            # Create the corresponding award.
            award = Award.objects.create(user=user, badge=badge)
            logger.info("award %s created for %s" % (award.badge.name, user.email))
