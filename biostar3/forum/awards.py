"""
Award specifications
"""
import logging
from .models import *
from functools import *
from itertools import *

logger = logging.getLogger('biostar')


class AwardDef(object):
    def __init__(self, uuid, name, desc, func, icon, type=Badge.BRONZE):
        self.uuid = uuid
        self.name = name
        self.desc = desc
        self.fun = func
        self.icon = icon
        self.template = "badge/default.html"
        self.type = type

        self.set_badge()

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

    def set_badge(self):
        self.badge = Badge.objects.filter(uuid=uuid).first()
        if not self.badge:
            self.badge = Badge.objects.create(uuid=self.uuid, name=self.name,
                                             desc=self.desc, icon=self.icon, type=self.type)


def wrap_list(obj, cond, func=lambda x: x):
    return func(obj) if cond else func()


def init_awards():

    #
    # Award definitions.
    # Do not change the uuids below after initializing the site! Doing so will create a new badge.
    #
    AUTOBIO = AwardDef(
        uuid='autobio',
        name="Autobiographer",
        desc="has more than 80 characters in the information field of the user's profile",
        func=lambda user: wrap_list(user, cond=len(user.profile.info) > 80),
        icon='<i class="fa fa-bullhorn"></i>'
    )

    GOOD_QUESTION = AwardDef(
        uuid = "goodquestion",
        name="Good Question",
        desc="asked a question that was upvoted at least 5 times",
        func=lambda user: Post.objects.filter(vote_count__gt=5, author=user, type=Post.QUESTION),
        icon="fa fa-question"
    )

    GOOD_ANSWER = AwardDef(
        uuid="goodanswer",
        name="Good Answer",
        desc="created an answer that was upvoted at least 5 times",
        func=lambda user: Post.objects.filter(vote_count__gt=5, author=user, type=Post.ANSWER),
        icon="fa fa-pencil-square-o"
    )

    STUDENT = AwardDef(
        uuid="student",
        name="Student",
        desc="asked a question with at least 3 up-votes",
        func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.QUESTION),
        icon="fa fa-certificate"
    )

    TEACHER = AwardDef(
        uuid='teacher',
        name="Teacher",
        desc="created an answer with at least 3 up-votes",
        func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.ANSWER),
        icon="fa fa-smile-o"
    )

    COMMENTATOR = AwardDef(
        uuid="commentator",
        name="Commentator",
        desc="created a comment with at least 3 up-votes",
        func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.COMMENT),
        icon="fa fa-comment"
    )

    CENTURION = AwardDef(
        uuid="centurion",
        name="Centurion",
        desc="created 100 posts",
        func=lambda user: wrap_list(user, Post.objects.filter(author=user).count() > 100),
        icon="fa fa-bolt",
        type=Badge.SILVER,
    )

    EPIC_QUESTION = AwardDef(
        uuid="epic-question",
        name="Epic Question",
        desc="created a question with more than 10,000 views",
        func=lambda user: Post.objects.filter(author=user, view_count__gt=10000),
        icon="fa fa-bullseye",
        type=Badge.GOLD,
    )

    POPULAR = AwardDef(
        uuid="popular-question",
        name="Popular Question",
        desc="created a question with more than 1,000 views",
        func=lambda user: Post.objects.filter(author=user, view_count__gt=1000),
        icon="fa fa-eye",
        type=Badge.GOLD,
    )

    ORACLE = AwardDef(
        uuid="oracle",
        name="Oracle",
        desc="created more than 1,000 posts (questions + answers + comments)",
        func=lambda user: wrap_list(user, Post.objects.filter(author=user).count() > 1000),
        icon="fa fa-sun-o",
        type=Badge.GOLD,
    )

    PUNDIT = AwardDef(
        uuid="pundit",
        name="Pundit",
        desc="created a comment with more than 10 votes",
        func=lambda user: Post.objects.filter(author=user, type=Post.COMMENT, vote_count__gt=10),
        icon="fa fa-comments-o",
        type=Badge.SILVER,
    )

    GURU = AwardDef(
        uuid="guru",
        name="Guru",
        desc="received more than 100 upvotes",
        func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 100),
        icon="fa fa-beer",
        type=Badge.SILVER,
    )

    CYLON = AwardDef(
        uuid="cylon",
        name="Cylon",
        desc="received 1,000 up votes",
        func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 1000),
        icon="fa fa-rocket",
        type=Badge.GOLD,
    )

    VOTER = AwardDef(
        uuid="voter",
        name="Voter",
        desc="voted more than 100 times",
        func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 100),
        icon="fa fa-thumbs-o-up"
    )

    SUPPORTER = AwardDef(
        uuid="supporter",
        name="Supporter",
        desc="voted at least 25 times",
        func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 25),
        icon="fa fa-thumbs-up",
        type=Badge.SILVER,
    )

    SCHOLAR = AwardDef(
        uuid="scholar",
        name="Scholar",
        desc="created an answer that has been accepted",
        func=lambda user: Post.objects.filter(author=user, type=Post.ANSWER, has_accepted=True),
        icon="fa fa-check-circle-o"
    )

    PROPHET = AwardDef(
        uuid="prophet",
        name="Prophet",
        desc="created a post with more than 20 followers",
        func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, subs_count__gt=20),
        icon="fa fa-pagelines"
    )

    LIBRARIAN = AwardDef(
        uuid="librarian",
        name="Librarian",
        desc="created a post with more than 10 bookmarks",
        func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, book_count__gt=10),
        icon="fa fa-bookmark-o"
    )

    awards = [
        AUTOBIO,
    ]


    return awards
