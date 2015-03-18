"""
Award specifications
"""
import logging

from django.db.models import Q, F

from .models import *
from . import mailer

from functools import *
from itertools import *

logger = logging.getLogger('biostar')

default_select_func = lambda: []
default_cond_func = lambda obj: True
default_date_func = lambda obj: right_now()


def find_vote_date(post, N=1):
    """
    Finds the date of the nth vote by date.
    """
    vote = Vote.objects.all().order_by("date")[N:N + 1].first()
    if not vote:
        return right_now()
    else:
        return vote.date


class AwardDef(object):
    cache = dict()

    def __init__(self, uuid, name, desc, icon,
                 type=Badge.BRONZE,
                 selector=default_select_func,
                 cond_func=default_cond_func,
                 date_func=default_date_func):

        self.uuid = uuid
        self.name = name
        self.desc = desc
        self.icon = icon
        self.type = type

        self.template = "badge/default.html"

        # The function that selects candidates for the badge.
        self.selector = selector

        # The function that determines the condition to award the badge
        self.cond_func = cond_func

        # The function that determines the date of the badge
        self.date_func = date_func

        # Initialize the corresponding badge.
        self.init_badge()

    def check(self, user, post=None, cache={}, override=False):
        """
        Attempts to create all awards based on the selector.
        """
        try:

            # Find qualifying targets
            targets = ifilter(self.cond_func, self.selector)

            for obj in targets:

                date = self.date_func(obj)

                if isinstance(obj, User):
                    # Object is a user.
                    user, post = obj, None
                else:
                    # Object must be a post.
                    user, post = obj.author, obj

                # Create the award
                award = Award.objects.create(badge=self.badge, user=user, post=post, date=date)

                # Generate the message for the award.
                em = mailer.EmailTemplate("award_created_message.html")

                # Create a local message
                em.create_messages(author=user, targets=[award])

                # TODO: Send an email message if the user is tracking group via email.
                return award

        except KeyError, exc:
            logger.error("award %s error %s" % (self.uuid, exc))

        return None

    def init_badge(self):
        """
        Initializes the badge underlying the award.
        """
        self.badge = Badge.objects.filter(uuid=self.uuid).first()
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
    def autobio_selector(uuid):
        cond = Q(profile__info='') | Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).select_related("profile")
        return query

    AUTOBIO = AwardDef(
        uuid='autobio',
        name="Autobiographer",
        desc="has more than 80 characters in the information field of the user's profile",
        selector=autobio_selector(uuid),
        cond_func=lambda user: len(user.profile.info) > 80,
        icon='<i class="fa fa-bullhorn"></i>'
    )

    def good_question_selector(uuid):
        cond = Q(vote_count__lt=5) | Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).select_related("author")
        return query

    GOOD_QUESTION = AwardDef(
        uuid="goodquestion",
        name="Good Question",
        desc="asked a question that was upvoted at least 5 times",
        selector = good_question_selector(uuid),
        date_func=lambda post: find_vote_date(post, 5),
        icon='<i class="fa fa-question"></i>',
    )

    '''
    GOOD_ANSWER = AwardDef(
        uuid="goodanswer",
        name="Good Answer",
        desc="created an answer that was upvoted at least 5 times",
        cond_func=lambda user, post: Post.objects.filter(vote_count__gt=5, author=user, type=Post.ANSWER),
        date_func=lambda user, post: find_vote_date(post, 5),
        icon="fa fa-pencil-square-o"
    )

    STUDENT = AwardDef(
        uuid="student",
        name="Student",
        desc="asked a question with at least 3 up-votes",
        # func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.QUESTION),
        icon="fa fa-certificate"
    )

    TEACHER = AwardDef(
        uuid='teacher',
        name="Teacher",
        desc="created an answer with at least 3 up-votes",
        # func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.ANSWER),
        icon="fa fa-smile-o"
    )

    COMMENTATOR = AwardDef(
        uuid="commentator",
        name="Commentator",
        desc="created a comment with at least 3 up-votes",
        # func=lambda user: Post.objects.filter(vote_count__gt=2, author=user, type=Post.COMMENT),
        icon="fa fa-comment"
    )

    CENTURION = AwardDef(
        uuid="centurion",
        name="Centurion",
        desc="created 100 posts",
        # func=lambda user: wrap_list(user, Post.objects.filter(author=user).count() > 100),
        icon="fa fa-bolt",
        type=Badge.SILVER,
    )

    EPIC_QUESTION = AwardDef(
        uuid="epic-question",
        name="Epic Question",
        desc="created a question with more than 10,000 views",
        # func=lambda user: Post.objects.filter(author=user, view_count__gt=10000),
        icon="fa fa-bullseye",
        type=Badge.GOLD,
    )

    POPULAR = AwardDef(
        uuid="popular-question",
        name="Popular Question",
        desc="created a question with more than 1,000 views",
        # func=lambda user: Post.objects.filter(author=user, view_count__gt=1000),
        icon="fa fa-eye",
        type=Badge.GOLD,
    )

    ORACLE = AwardDef(
        uuid="oracle",
        name="Oracle",
        desc="created more than 1,000 posts (questions + answers + comments)",
        # func=lambda user: wrap_list(user, Post.objects.filter(author=user).count() > 1000),
        icon="fa fa-sun-o",
        type=Badge.GOLD,
    )

    PUNDIT = AwardDef(
        uuid="pundit",
        name="Pundit",
        desc="created a comment with more than 10 votes",
        # func=lambda user: Post.objects.filter(author=user, type=Post.COMMENT, vote_count__gt=10),
        icon="fa fa-comments-o",
        type=Badge.SILVER,
    )

    GURU = AwardDef(
        uuid="guru",
        name="Guru",
        desc="received more than 100 upvotes",
        # func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 100),
        icon="fa fa-beer",
        type=Badge.SILVER,
    )

    CYLON = AwardDef(
        uuid="cylon",
        name="Cylon",
        desc="received 1,000 up votes",
        # func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 1000),
        icon="fa fa-rocket",
        type=Badge.GOLD,
    )

    VOTER = AwardDef(
        uuid="voter",
        name="Voter",
        desc="voted more than 100 times",
        # func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 100),
        icon="fa fa-thumbs-o-up"
    )

    SUPPORTER = AwardDef(
        uuid="supporter",
        name="Supporter",
        desc="voted at least 25 times",
        # func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 25),
        icon="fa fa-thumbs-up",
        type=Badge.SILVER,
    )

    SCHOLAR = AwardDef(
        uuid="scholar",
        name="Scholar",
        desc="created an answer that has been accepted",
        # func=lambda user: Post.objects.filter(author=user, type=Post.ANSWER, has_accepted=True),
        icon="fa fa-check-circle-o"
    )

    PROPHET = AwardDef(
        uuid="prophet",
        name="Prophet",
        desc="created a post with more than 20 followers",
        # func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, subs_count__gt=20),
        icon="fa fa-pagelines"
    )

    LIBRARIAN = AwardDef(
        uuid="librarian",
        name="Librarian",
        desc="created a post with more than 10 bookmarks",
        # func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, book_count__gt=10),
        icon="fa fa-bookmark-o"
    )
    '''
    awards = [
        AUTOBIO,
        GOOD_QUESTION,
    ]

    return awards
