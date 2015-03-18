"""
Award specifications
"""
import logging

from django.db.models import Q, F, Count, Avg

from .models import *
from . import mailer

from functools import *
from itertools import *

logger = logging.getLogger('biostar')

default_selector = lambda: []
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

# How many posts to qualify for centurion.
CENTURION_COUNT = 100


class AwardDef(object):
    cache = dict()

    def __init__(self, uuid, name, desc, icon,
                 type=Badge.BRONZE,
                 selector=default_selector,
                 date_func=default_date_func):

        self.uuid = uuid
        self.name = name
        self.desc = desc
        self.icon = icon
        self.type = type

        self.template = "badge/default.html"

        # The function that selects candidates for the badge.
        self.selector = selector

        # The function that determines the date of the badge
        self.date_func = date_func

        # Initialize the corresponding badge.
        self.init_badge()

    def check(self, user, post=None, cache={}, override=False):
        """
        Attempts to create all awards based on the selector.
        """
        try:

            targets = self.selector(self.uuid)

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

                print "creating award %s" % award.badge.uuid

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
        query = User.objects.exclude(cond).select_related("profile").distinct()
        query = ifilter(lambda user: len(user.profile.info) > 80, query)
        return query

    AUTOBIO = AwardDef(
        uuid='autobio',
        name="Autobiographer",
        desc="has more than 80 characters in the information field of the user's profile",
        selector=autobio_selector,
        icon='<i class="fa fa-bullhorn"></i>'
    )

    def good_question_selector(uuid):
        cond = Q(vote_count__lt=5) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.QUESTION).exclude(cond).select_related("author").distinct()
        return query

    GOOD_QUESTION = AwardDef(
        uuid="goodquestion",
        name="Good Question",
        desc="asked a question that was upvoted at least 5 times",
        selector=good_question_selector,
        date_func=lambda post: find_vote_date(post, 5),
        icon='<i class="fa fa-question"></i>',
    )

    def good_answer_selector(uuid):
        cond = Q(vote_count__lt=5) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.ANSWER).exclude(cond).select_related("author")
        return query

    GOOD_ANSWER = AwardDef(
        uuid="goodanswer",
        name="Good Answer",
        desc="created an answer that was upvoted at least 5 times",
        selector=good_answer_selector,
        date_func=lambda post: find_vote_date(post, 5),
        icon='<i class="fa fa-pencil-square-o"></i>'
    )

    def student_selector(uuid):
        cond = Q(vote_count__lt=3) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.QUESTION).exclude(cond).distinct()
        return query

    STUDENT = AwardDef(
        uuid="student",
        name="Student",
        desc="asked a question with at least 3 up-votes",
        selector=student_selector,
        icon='<i class="fa fa-certificate"></i>'
    )

    def teacher_selector(uuid):
        cond = Q(vote_count__lt=3) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.ANSWER).exclude(cond).distinct()
        return query

    TEACHER = AwardDef(
        uuid='teacher',
        name="Teacher",
        desc="created an answer with at least 3 up-votes",
        selector=teacher_selector,
        icon='<i class="fa fa-smile-o"></i>'
    )

    def commentator_selector(uuid):
        cond = Q(vote_count__lt=3) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.COMMENT).exclude(cond).distinct()
        return query


    COMMENTATOR = AwardDef(
        uuid="commentator",
        name="Commentator",
        desc="created a comment with at least 3 up-votes",
        selector=commentator_selector,
        icon='<i class="fa fa-comment-o"></i>'
    )

    def centurion_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).annotate(count=Count('post')).filter(count__gte=CENTURION_COUNT).distinct()
        return query

    CENTURION = AwardDef(
        uuid="centurion",
        name="Centurion",
        desc="created more than 100 posts",
        selector=centurion_selector,
        icon='<i class="fa fa-bolt"></i>',
        type=Badge.SILVER,
    )

    '''
    EPIC_QUESTION = AwardDef(
        uuid="epic-question",
        name="Epic Question",
        desc="created a question with more than 10,000 views",
        # func=lambda user: Post.objects.filter(author=user, view_count__gt=10000),
        icon='<i class="fa fa-bullseye"></i>'
        type=Badge.GOLD,
    )

    POPULAR = AwardDef(
        uuid="popular-question",
        name="Popular Question",
        desc="created a question with more than 1,000 views",
        # func=lambda user: Post.objects.filter(author=user, view_count__gt=1000),
        icon='<i class="fa fa-eye"></i>'
        type=Badge.GOLD,
    )

    ORACLE = AwardDef(
        uuid="oracle",
        name="Oracle",
        desc="created more than 1,000 posts (questions + answers + comments)",
        # func=lambda user: wrap_list(user, Post.objects.filter(author=user).count() > 1000),
        icon='<i class="fa fa-sun-o"></i>'
        type=Badge.GOLD,
    )

    PUNDIT = AwardDef(
        uuid="pundit",
        name="Pundit",
        desc="created a comment with more than 10 votes",
        # func=lambda user: Post.objects.filter(author=user, type=Post.COMMENT, vote_count__gt=10),
        icon='<i class="fa fa-comments-o"></i>'
        type=Badge.SILVER,
    )

    GURU = AwardDef(
        uuid="guru",
        name="Guru",
        desc="received more than 100 upvotes",
        # func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 100),
        icon='<i class="fa fa-beer"></i>'
        type=Badge.SILVER,
    )

    CYLON = AwardDef(
        uuid="cylon",
        name="Cylon",
        desc="received 1,000 up votes",
        # func=lambda user: wrap_list(user, Vote.objects.filter(post__author=user).count() > 1000),
        icon='<i class="fa fa-rocket"></i>',
        type=Badge.GOLD,
    )

    VOTER = AwardDef(
        uuid="voter",
        name="Voter",
        desc="voted more than 100 times",
        # func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 100),
        icon='<i class="fa fa-thumbs-o"></i>',
    )

    SUPPORTER = AwardDef(
        uuid="supporter",
        name="Supporter",
        desc="voted at least 25 times",
        # func=lambda user: wrap_list(user, Vote.objects.filter(author=user).count() > 25),
        icon='<i class="fa fa-thumbs-up"></i>',
        type=Badge.SILVER,
    )

    SCHOLAR = AwardDef(
        uuid="scholar",
        name="Scholar",
        desc="created an answer that has been accepted",
        # func=lambda user: Post.objects.filter(author=user, type=Post.ANSWER, has_accepted=True),
        icon='<i class="fa fa-check-circle-o"></i>',
    )

    PROPHET = AwardDef(
        uuid="prophet",
        name="Prophet",
        desc="created a post with more than 20 followers",
        # func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, subs_count__gt=20),
        icon='<i class="fa fa-pagelines"></i>',
    )

    LIBRARIAN = AwardDef(
        uuid="librarian",
        name="Librarian",
        desc="created a post with more than 10 bookmarks",
        # func=lambda user: Post.objects.filter(author=user, type__in=Post.TOP_LEVEL, book_count__gt=10),
        icon='<i class="fa fa-bookmark-o"></i>',
    )
    '''
    awards = [
        AUTOBIO,
        GOOD_QUESTION,
        GOOD_ANSWER,
        STUDENT,
        TEACHER,
        COMMENTATOR,
        CENTURION,
    ]

    return awards
