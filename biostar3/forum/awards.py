"""
Award specifications
"""

from django.db.models import Q, Count

from .models import *
from . import mailer
from biostar3.utils.compat import *

logger = logging.getLogger('biostar')

default_selector = lambda: []
default_date_func = lambda obj: right_now()
default_user_date_func = lambda obj: obj.profile.date_joined
default_post_date_func = lambda obj: obj.creation_date

# How many posts to qualify for centurion.
CENTURION_COUNT = 100
EPIC_COUNT = 10000
POPULAR_COUNT = 1000
ORACLE_COUNT = 1000
PUNDIT_COUNT = 10
GURU_COUNT = 100
CYLON_COUNT = 1000
VOTER_COUNT = 100
SUPPORTER_COUNT = 25
PROPHET_COUNT = 20
BOOKMARK_COUNT = 20


class AwardDef(object):
    cache = dict()

    def __init__(self, uuid, name, desc, icon,
                 type=Badge.BRONZE,
                 selector=default_selector,
                 date_func=None):

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

    def check(self, cache={}):
        """
        Attempts to create all awards based on the selector.
        """
        try:

            targets = self.selector(self.uuid)

            for obj in targets:

                if isinstance(obj, User):
                    # Object is a user.
                    user, post = obj, None
                    date = default_user_date_func(user)
                else:
                    # Object must be a post.
                    user, post = obj.author, obj
                    date = default_post_date_func(post)

                # Overrride the date function if necessary.
                if self.date_func:
                    date = self.date_func(obj)

                # Create the award
                award = Award.objects.create(badge=self.badge, user=user, post=post, date=date)

                # Generate the message for the award.
                data = dict(award=award, badge=award.badge, post=post, user=user)
                em = mailer.EmailTemplate("award_created_message.html", data=data)

                # Create a local message
                em.create_messages(author=user, users=[award.user])

            self.badge.count = Award.objects.filter(badge=self.badge).count()
            self.badge.save()

        except KeyError as exc:
            logger.error("award %s error %s" % (self.uuid, exc))


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


def get_awards():
    #
    # Award definitions.
    # Do not change the uuids below after initializing the site! Doing so will create a new badge.
    #
    def autobio_selector(uuid):
        cond = Q(profile__info='') | Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).select_related("profile").distinct()
        query = filter(lambda user: len(user.profile.info) > 80, query)
        return query

    AUTOBIO = AwardDef(
        uuid='autobio',
        name="Autobiographer",
        desc="Has more than 80 characters in the information field of the user's profile.",
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
        desc="Asked a question that was upvoted at least 5 times.",
        selector=good_question_selector,
        icon='<i class="fa fa-question"></i>',
    )

    def good_answer_selector(uuid):
        cond = Q(vote_count__lt=5) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.ANSWER).exclude(cond).select_related("author")
        return query

    GOOD_ANSWER = AwardDef(
        uuid="goodanswer",
        name="Good Answer",
        desc="Created an answer that was upvoted at least 5 times.",
        selector=good_answer_selector,
        icon='<i class="fa fa-pencil-square-o"></i>'
    )

    def student_selector(uuid):
        cond = Q(vote_count__lt=3) | Q(award__badge__uuid=uuid)
        query = Post.objects.filter(type=Post.QUESTION).exclude(cond).distinct()
        return query

    STUDENT = AwardDef(
        uuid="student",
        name="Student",
        desc="Asked a question with at least 3 up-votes.",
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
        desc="Created an answer with at least 3 up-votes.",
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
        desc="Created a comment with at least 3 up-votes.",
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
        desc="Created more than 100 posts.",
        selector=centurion_selector,
        icon='<i class="fa fa-bolt"></i>',
        type=Badge.SILVER,
    )

    def epic_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(type=Post.QUESTION, view_count__gte=EPIC_COUNT).distinct()
        return query

    EPIC_QUESTION = AwardDef(
        uuid="epic-question",
        name="Epic Question",
        desc="Created a question with more than 10,000 views.",
        selector=epic_selector,
        icon='<i class="fa fa-bullseye"></i>',
        type=Badge.GOLD,
    )

    def popular_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(type=Post.QUESTION, view_count__gte=POPULAR_COUNT).distinct()
        return query

    POPULAR_QUESTION = AwardDef(
        uuid="popular-question",
        name="Popular Question",
        desc="Created a question with more than 1,000 views.",
        selector=popular_selector,
        icon='<i class="fa fa-eye"></i>',
        type=Badge.GOLD,
    )

    def oracle_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).annotate(count=Count('post')).filter(count__gte=ORACLE_COUNT).distinct()
        return query

    ORACLE = AwardDef(
        uuid="oracle",
        name="Oracle",
        desc="Created more than 1,000 posts (questions + answers + comments).",
        selector=oracle_selector,
        icon='<i class="fa fa-sun-o"></i>',
        type=Badge.GOLD,
    )

    def pundit_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(vote_count__gte=PUNDIT_COUNT, type=Post.COMMENT).distinct()
        return query

    PUNDIT = AwardDef(
        uuid="pundit",
        name="Pundit",
        desc="Created a comment with more than 10 votes.",
        selector=pundit_selector,
        icon='<i class="fa fa-comments-o"></i>',
        type=Badge.SILVER,
    )

    def guru_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).filter(score__gte=GURU_COUNT).distinct()
        return query

    GURU = AwardDef(
        uuid="guru",
        name="Guru",
        desc="Received more than 100 upvotes.",
        selector=guru_selector,
        icon='<i class="fa fa-beer"></i>',
        type=Badge.SILVER,
    )

    def cylon_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).filter(score__gte=CYLON_COUNT).distinct()
        return query

    CYLON = AwardDef(
        uuid="cylon",
        name="Cylon",
        desc="Received more than 1,000 up votes.",
        selector=cylon_selector,
        icon='<i class="fa fa-rocket"></i>',
        type=Badge.GOLD,
    )

    def voter_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).annotate(count=Count("vote")).filter(count__gte=VOTER_COUNT).distinct()
        return query

    VOTER = AwardDef(
        uuid="voter",
        name="Voter",
        desc="Voted more than 100 times.",
        selector=voter_selector,
        icon='<i class="fa fa-thumbs-o"></i>',
    )

    def supporter_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = User.objects.exclude(cond).annotate(count=Count("vote")).filter(count__gte=SUPPORTER_COUNT).distinct()
        return query

    SUPPORTER = AwardDef(
        uuid="supporter",
        name="Supporter",
        desc="Voted at least 25 times.",
        selector=supporter_selector,
        icon='<i class="fa fa-thumbs-up"></i>',
        type=Badge.SILVER,
    )


    def scholar_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(has_accepted=True, type=Post.ANSWER).distinct()
        return query

    SCHOLAR = AwardDef(
        uuid="scholar",
        name="Scholar",
        desc="Created an answer that has been accepted.",
        selector=scholar_selector,
        icon='<i class="fa fa-check-circle-o"></i>',
    )

    def prophet_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(type__in=Post.TOP_LEVEL, subs_count__gte=PROPHET_COUNT).distinct()
        return query

    PROPHET = AwardDef(
        uuid="prophet",
        name="Prophet",
        desc="Created a post with more than 20 followers.",
        selector=prophet_selector,
        icon='<i class="fa fa-pagelines"></i>',
    )

    def librarian_selector(uuid):
        cond = Q(award__badge__uuid=uuid)
        query = Post.objects.exclude(cond).filter(book_count__gte=BOOKMARK_COUNT).distinct()
        return query

    LIBRARIAN = AwardDef(
        uuid="librarian",
        name="Librarian",
        desc="Created a post with more than 10 bookmarks.",
        selector=librarian_selector,
        icon='<i class="fa fa-bookmark-o"></i>',
    )

    awards = [
        AUTOBIO,
        GOOD_QUESTION,
        GOOD_ANSWER,
        STUDENT,
        TEACHER,
        COMMENTATOR,
        CENTURION,
        EPIC_QUESTION,
        POPULAR_QUESTION,
        ORACLE,
        PUNDIT,
        GURU,
        CYLON,
        VOTER,
        SUPPORTER,
        SCHOLAR,
        PROPHET,
        LIBRARIAN,
    ]

    return awards
