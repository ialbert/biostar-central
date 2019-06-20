import logging

from django.utils.timezone import utc
from datetime import datetime, timedelta

from biostar.accounts.models import User
from biostar.forum.models import Post, Vote, Badge, Award

logger = logging.getLogger("engine")


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def wrap_qs(cond, klass, pk):

    return klass.objects.filter(pk=pk) if cond else klass.objects.none()


class AwardDef(object):
    def __init__(self, name, desc, func, icon, max_awarded=None, type=Badge.BRONZE):
        self.name = name
        self.desc = desc
        self.fun = func
        self.icon = icon
        self.template = ""
        self.type = type
        # Max number of times this award can be given
        # No limit if left empty.
        self.max_awarded = max_awarded

    def validate(self, *args, **kwargs):

        user = args[0]
        try:
            value = self.fun(*args, **kwargs).order_by("pk")
        except Exception as exc:
            logger.error("validator error %s" % exc)
            return []

        # Award user has won at this point.
        awarded = Award.objects.filter(badge__name=self.name, user=user)

        # Exclude targets that already awarded.
        if isinstance(value, Post):
            value = value.exclude(id__in=awarded.post_set)
        else:
            value = value.exclude(pk=user.id)

        # Ensure users does not get over awarded.
        if self.max_awarded and len(awarded) >= self.max_awarded:
            return []

        return value

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name


# Award definitions
AUTOBIO = AwardDef(
    name="Autobiographer",
    desc="has more than 80 characters in the information field of the user's profile",
    func=lambda user: wrap_qs(len(user.profile.text) > 80, User, user.id),
    max_awarded=1,
    icon="bullhorn icon"
)

GOOD_QUESTION = AwardDef(
    name="Good Question",
    desc="asked a question that was upvoted at least 5 times",
    func=lambda user: Post.objects.filter(vote_count__gte=5, author=user, type=Post.QUESTION),
    max_awarded=1,
    icon="question icon"
)

GOOD_ANSWER = AwardDef(
    name="Good Answer",
    desc="created an answer that was upvoted at least 5 times",
    func=lambda user: Post.objects.filter(vote_count__gt=5, author=user, type=Post.ANSWER),
    max_awarded=1,
    icon="edit outline icon"
)

STUDENT = AwardDef(
    name="Student",
    desc="asked a question with at least 3 up-votes",
    func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.QUESTION),
    max_awarded=1,
    icon="certificate icon"
)

TEACHER = AwardDef(
    name="Teacher",
    desc="created an answer with at least 3 up-votes",
    func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.ANSWER),
    max_awarded=1,
    icon="smile outline icon"
)

COMMENTATOR = AwardDef(
    name="Commentator",
    desc="created a comment with at least 3 up-votes",
    func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.COMMENT),
    max_awarded=1,
    icon="comment icon"
)

CENTURION = AwardDef(
    name="Centurion",
    desc="created 100 posts",
    func=lambda user: wrap_qs(Post.objects.filter(author=user).count() > 100, User, user.id),
    max_awarded=1,
    icon="bolt icon",
    type=Badge.SILVER,
)

EPIC_QUESTION = AwardDef(
    name="Epic Question",
    desc="created a question with more than 10,000 views",
    func=lambda user: Post.objects.filter(author=user, view_count__gt=10000),
    max_awarded=1,
    icon="bullseye icon",
    type=Badge.GOLD,
)

POPULAR = AwardDef(
    name="Popular Question",
    desc="created a question with more than 1,000 views",
    func=lambda user: Post.objects.filter(author=user, view_count__gt=1000),
    max_awarded=1,
    icon="eye icon",
    type=Badge.GOLD,
)

ORACLE = AwardDef(
    name="Oracle",
    desc="created more than 1,000 posts (questions + answers + comments)",
    func=lambda user: wrap_qs(Post.objects.filter(author=user).count() > 1000, User, user.id),
    max_awarded=1,
    icon="sun icon",
    type=Badge.GOLD,
)

PUNDIT = AwardDef(
    name="Pundit",
    desc="created a comment with more than 10 votes",
    func=lambda user: Post.objects.filter(author=user, type=Post.COMMENT, vote_count__gt=10),
    max_awarded=1,
    icon="comments icon",
    type=Badge.SILVER,
)

GURU = AwardDef(
    name="Guru",
    desc="received more than 100 upvotes",
    func=lambda user: wrap_qs(Vote.objects.filter(post__author=user).count() > 100, User, user.id),
    max_awarded=1,
    icon="beer icon",
    type=Badge.SILVER,
)

CYLON = AwardDef(
    name="Cylon",
    desc="received 1,000 up votes",
    func=lambda user: wrap_qs(Vote.objects.filter(post__author=user).count() > 1000, User, user.id),
    max_awarded=1,
    icon="rocket icon",
    type=Badge.GOLD,
)

VOTER = AwardDef(
    name="Voter",
    desc="voted more than 100 times",
    func=lambda user: wrap_qs(Vote.objects.filter(author=user).count() > 100, User, user.id),
    max_awarded=1,
    icon="thumbs up outline"
)

SUPPORTER = AwardDef(
    name="Supporter",
    desc="voted at least 25 times",
    func=lambda user: wrap_qs(Vote.objects.filter(author=user).count() > 25, User, user.id),
    max_awarded=1,
    icon="thumbs up icon",
    type=Badge.SILVER,
)

SCHOLAR = AwardDef(
    name="Scholar",
    desc="created an answer that has been accepted",
    func=lambda user: Post.objects.filter(author=user, type=Post.ANSWER, accept_count__gt=0),
    max_awarded=1,
    icon="check circle outline icon"
)

PROPHET = AwardDef(
    name="Prophet",
    desc="created a post with more than 20 followers",
    func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, subs_count__gt=20),
    max_awarded=1,
    icon="leaf icon"
)

LIBRARIAN = AwardDef(
    name="Librarian",
    desc="created a post with more than 10 bookmarks",
    func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, book_count__gt=10),
    max_awarded=1,
    icon="bookmark outline icon"
)


def rising_star(user):
    # The user joined no more than three months ago
    cond = now() < user.profile.date_joined + timedelta(weeks=15)
    cond = cond and Post.objects.filter(author=user).count() > 50
    return wrap_qs(cond, User, user.id)

RISING_STAR = AwardDef(
    name="Rising Star",
    desc="created 50 posts within first three months of joining",
    func=rising_star,
    icon="star icon",
    max_awarded=1,
    type=Badge.GOLD,
)


GREAT_QUESTION = AwardDef(
    name="Great Question",
    desc="created a question with more than 5,000 views",
    func=lambda user: Post.objects.filter(author=user, view_count__gt=5000),
    icon="fire icon",
    type=Badge.SILVER,
)

GOLD_STANDARD = AwardDef(
    name="Gold Standard",
    desc="created a post with more than 25 bookmarks",
    func=lambda user: Post.objects.filter(author=user, book_count__gt=25),
    icon="bookmark icon",
    type=Badge.GOLD,
)

APPRECIATED = AwardDef(
    name="Appreciated",
    desc="created a post with more than 5 votes",
    func=lambda user: Post.objects.filter(author=user, vote_count__gt=4),
    icon="heart icon",
    type=Badge.SILVER,
)


ALL_AWARDS = [

    # These awards can only be earned once
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
    EPIC_QUESTION,
    ORACLE,
    PUNDIT,
    GOOD_ANSWER,
    GOOD_QUESTION,
    PROPHET,
    LIBRARIAN,
    # These awards can be won multiple times
    GREAT_QUESTION,
    GOLD_STANDARD,
    APPRECIATED,
]
