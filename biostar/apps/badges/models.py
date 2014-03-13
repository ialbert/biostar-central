from django.db import models
from django.conf import settings
from django.core.urlresolvers import reverse
from django.core import mail

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

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)


# Simple badges can be computed via a single function call
AUTOBIO = BadgeDef(
    name = "Autobiographer",
    desc = "more than 80 character in the information field of your profile",
    func = lambda user: (len(user.profile.info) > 80)
)

USER_BADGES = [
    AUTOBIO,
]


# Tries to award a badge to the user
def check_badges(request):
    user = request.user
    awards = set(a.name for a in Award.objects.filter(user=user))

    for obj in USER_BADGES:
        if obj.fun(user) and obj.name not in awards:
            badge = Badge.objects.get_or_create(name=obj.name)
            award = Award.objects.create(user=user, badge=badge)
            award.notify_user()
            return

