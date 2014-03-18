from django.db import models
from django.conf import settings
from django.core.urlresolvers import reverse
from django.core import mail
import logging

from django.utils.timezone import utc
from datetime import datetime, timedelta

def now():
    return datetime.utcnow().replace(tzinfo=utc)


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
    date = models.DateTimeField(auto_now_add=True)

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
            return int(value)
        except Exception, exc:
            logger.error("validator error %s" % exc)
        return 0

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)


from biostar.apps.posts.models import Post, Vote

# Award definitions
AUTOBIO = AwardDef(
    name = "Autobiographer",
    desc = "has more than 80 characters in the information field of the user's profile",
    func = lambda user: (len(user.profile.info) > 80),
    icon = "fa fa-bullhorn"
)

GOOD_QUESTION = AwardDef(
    name = "Good Question",
    desc = "asked a question that was upvoted",
    func = lambda user: Post.objects.filter(vote_count__gt=0, author=user, type=Post.QUESTION).count() > 0,
    icon = "fa fa-question"
)

GOOD_ANSWER = AwardDef(
    name = "Good Answer",
    desc = "created an answer that was upvoted",
    func = lambda user: Post.objects.filter(vote_count__gt=0, author=user, type=Post.ANSWER).count() > 0,
    icon = "fa fa-pencil-square-o"
)

STUDENT = AwardDef(
    name = "Student",
    desc = "asked a question with at least three up-votes",
    func = lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.QUESTION).count(),
    icon = "fa fa-certificate"
)

TEACHER = AwardDef(
    name = "Teacher",
    desc = "created an answer with at least three up-votes",
    func = lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.ANSWER).count(),
    icon = "fa fa-smile-o"
)

COMMENTATOR = AwardDef(
    name = "Commentator",
    desc = "created a comment with at least one up-vote",
    func = lambda user: Post.objects.filter(vote_count__gt=0, author=user, type=Post.COMMENT).count(),
    icon = "fa fa-comment"
)

CENTURION = AwardDef(
    name = "Centurion",
    desc = "created 100 posts",
    func = lambda user: Post.objects.filter(author=user).count()>100,
    icon = "fa fa-bolt",
    type = Badge.SILVER,
)

BULLSEYE = AwardDef(
    name = "Bullseye",
    desc = "created a question with more than 10,000 views",
    func = lambda user: Post.objects.filter(author=user, view_count__gt=10000).count(),
    icon = "fa fa-bullseye",
    type = Badge.GOLD,
)

POPULAR = AwardDef(
    name = "Popular Question",
    desc = "created a question with more than 1,000 views",
    func = lambda user: Post.objects.filter(author=user, view_count__gt=1000).count(),
    icon = "fa fa-eye",
    type = Badge.GOLD,
)

SUPERNOVA = AwardDef(
    name = "Supernova",
    desc = "created more than 1,000 posts (questions + answers + comments)",
    func = lambda user: Post.objects.filter(author=user).count() > 1000,
    icon = "fa fa-sun-o",
    type = Badge.GOLD,
)

PUNDIT = AwardDef(
    name = "Pundit",
    desc = "created a comment with more than 10 votes",
    func = lambda user: Post.objects.filter(author=user, type=Post.COMMENT, vote_count__gt=10).count(),
    icon = "fa fa-comments-o",
    type = Badge.SILVER,
)

GURU = AwardDef(
    name = "Guru",
    desc = "received 50 upvotes",
    func = lambda user: Vote.objects.filter(post__author=user).count() > 50,
    icon = "fa fa-beer",
    type = Badge.SILVER,
)

CYLON = AwardDef(
    name = "Cylon",
    desc = "received 1,000 up votes",
    func = lambda user: Vote.objects.filter(post__author=user).count() > 1000,
    icon = "fa fa-rocket",
    type = Badge.GOLD,
)

VOTER = AwardDef(
    name = "Voter",
    desc = "upvoted one post",
    func = lambda user: Vote.objects.filter(author=user).count() > 0,
    icon = "fa fa-thumbs-o-up"
)

SUPPORTER = AwardDef(
    name = "Supporter",
    desc = "voted at least 25 times",
    func = lambda user: Vote.objects.filter(author=user).count() > 25,
    icon = "fa fa-thumbs-up",
    type = Badge.SILVER,
)

SCHOLAR = AwardDef(
    name = "Scholar",
    desc = "created an answer that has been accepted",
    func = lambda user: Post.objects.filter(author=user, type=Post.ANSWER, has_accepted=True).count() > 0,
    icon = "fa fa-check-circle-o"
)

def rising_star(user):
    # The user joined no more than three months ago
    cond = now() < user.profile.date_joined + timedelta(weeks=15)
    cond = cond and Post.objects.filter(author=user).count() > 50
    return cond

RISING_STAR = AwardDef(
    name = "Rising Star",
    desc = "created 50 posts within first three months of joining",
    func = rising_star,
    icon = "fa fa-star",
    type = Badge.GOLD,
)

# These awards can only be earned once
SINGLE_AWARDS = [
    AUTOBIO,
    STUDENT,
    TEACHER,
    COMMENTATOR,
    SUPPORTER,
    SCHOLAR,
    VOTER,
    CENTURION,
    CYLON,
    RISING_STAR,
    GURU,
    POPULAR,
    BULLSEYE,
    SUPERNOVA,
    PUNDIT,
    GOOD_ANSWER,
    GOOD_QUESTION,
]

def great_question(user):
    count = Post.objects.filter(author=user, view_count__gt=1000).count()
    return count

GREAT_QUESTION = AwardDef(
    name = "Great Question",
    desc = "created a question with 1,000 views",
    func = great_question,
    icon = "fa fa-fire",
    type = Badge.SILVER,
)

def oracle(user):
    votes = Vote.objects.filter(post__author=user, type=Vote.BOOKMARK)
    # find distinct posts, not supported on all backends
    posts = set( vote.post_id for vote in votes )
    return len(posts)

ORACLE = AwardDef(
    name = "Oracle",
    desc = "created a post with more than 25 bookmarks",
    func = oracle,
    icon = "fa fa-bookmark",
    type = Badge.GOLD,
)

APPRECIATED = AwardDef(
    name = "Appreciated",
    desc = "created a post with more than 5 votes",
    func = lambda user: Post.objects.filter(author=user, vote_count__gt=4).count(),
    icon = "fa fa-heart",
    type = Badge.SILVER,
)


# These awards can be won multiple times
MULTI_AWARDS = [
    GREAT_QUESTION,
    ORACLE,
    APPRECIATED,
]

ALL_AWARDS = SINGLE_AWARDS + MULTI_AWARDS