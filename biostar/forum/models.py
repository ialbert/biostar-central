import bleach
import logging
import datetime
import mistune

from django.utils import timezone
from django.db import models
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.sites.models import Site
from django.dispatch import receiver
from django.db.models.signals import post_save
from django.db.models import F, Q
from biostar.utils.shortcuts import reverse
from taggit.managers import TaggableManager
from biostar.engine.models import Project
from . import util

User = get_user_model()


# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN

logger = logging.getLogger("engine")


def get_sentinel_user():
    return User.objects.get_or_create(username='deleted').first()


class SubscriptionManager(models.Manager):
    def get_subs(self, post):
        "Returns all subscriptions for a post, exclude the "
        return self.filter(post=post.root).select_related("user")


class PostManager(models.Manager):

    def get_queryset(self):
        "Regular queries exclude deleted stuff"
        query = super().get_queryset().filter(project=None).prefetch_related("tags", "thread_users__profile",
                                                                             "thread_users")

        return query

    def get_discussions(self, **kwargs):
        "Get posts exclusively tied to projects"

        query = super().get_queryset().exclude(project=None, status=Post.DELETED).filter(**kwargs)
        query = query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile", "project")
        query = query.prefetch_related("tags", "thread_users")
        return query

    def get_all(self, **kwargs):
        "Return everything"
        query = super().get_queryset().filter(**kwargs)
        query = query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile")
        query = query.prefetch_related("tags", "thread_users", "thread_users__profile")
        return query

    def following(self, user):
        query = self.filter(~Q(subs__type=Subscription.NO_MESSAGES), subs__user=user).exclude(status=Post.DELETED)
        query = query.select_related("root", "author", "lastedit_user")
        query = query.prefetch_related("tags")
        return query

    def my_bookmarks(self, user):
        query = self.filter(votes__author=user, votes__type=Vote.BOOKMARK)
        query = query.select_related("root", "author", "lastedit_user")
        query = query.prefetch_related("tags")
        return query

    def my_post_votes(self, user):
        "Posts that received votes from other people "
        vote_query = Vote.objects.exclude(author=user).filter(post__in=self.filter(author=user))
        query = self.filter(author=user, votes__in=vote_query)
        query = query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile")
        query = query.prefetch_related("tags")
        return query

    def my_posts(self, target, user):

        if user.is_anonymous or target.is_anonymous:
            return self.filter().exclude(status=Post.DELETED)

        # Show all posts for moderators or targets
        if user.profile.is_moderator or user == target:
            query = self.filter(author=target)
        else:
            query = self.filter(author=target).exclude(status=Post.DELETED)

        query = query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile")
        query = query.prefetch_related("tags")
        query = query.order_by("-creation_date")
        return query

    def fixcase(self, text):
        return text.upper() if len(text) == 1 else text.lower()

    def tag_search(self, text, defer_content=True):
        "Performs a query by one or more , separated tags"
        include, exclude = [], []
        # Split the given tags on ',' and '+'.
        terms = text.split(',') if ',' in text else text.split('+')
        for term in terms:
            term = term.strip()
            if term.endswith("!"):
                exclude.append(self.fixcase(term[:-1]))
            else:
                include.append(self.fixcase(term))

        if include:
            query = self.filter(type__in=Post.TOP_LEVEL, tags__name__in=include).exclude(
                tags__name__in=exclude)
        else:
            query = self.filter(type__in=Post.TOP_LEVEL).exclude(tags__name__in=exclude)

        query = query.filter(status=Post.OPEN)
        if defer_content:
            # Remove fields that are not used.
            query = query.defer('content', 'html')

        # Get the tags.
        query = query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile").prefetch_related("tags").distinct()

        return query

    def get_thread(self, root, user):
        # Populate the object to build a tree that contains all posts in the thread.
        is_moderator = user.is_authenticated and user.profile.is_moderator
        if is_moderator:
            query = self.get_all(root=root).select_related("root", "parent", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile").order_by("type", "-has_accepted",
                                                                                                      "-vote_count", "creation_date")
        else:
            query = self.get_all(root=root).exclude(status=Post.DELETED).select_related("root", "parent", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile").order_by("type",
                                                                                         "-has_accepted",
                                                                                         "-vote_count",
                                                                                         "creation_date")


        return query

    def top_level(self, user):
        "Returns posts based on a user type"
        is_moderator = user.is_authenticated and (user.profile.is_moderator or user.profile.is_manager)
        if is_moderator:
            query = self.filter(type__in=Post.TOP_LEVEL)
        else:
            query = self.filter(type__in=Post.TOP_LEVEL).exclude(status=Post.DELETED)

        return query.select_related("root", "author", "author__profile",
                                    "lastedit_user", "lastedit_user__profile").prefetch_related("tags")#.defer("content", "html")


class Post(models.Model):
    "Represents a post in a forum"

    objects = PostManager()

    # Post statuses.
    PENDING, OPEN, CLOSED, DELETED = range(4)
    STATUS_CHOICES = [(PENDING, "Pending"), (OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types. Answers should be listed before comments.
    QUESTION, ANSWER, JOB, FORUM, PAGE, BLOG, COMMENT, DATA, TUTORIAL, BOARD, TOOL, NEWS = range(12)

    TYPE_CHOICES = [
        (QUESTION, "Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (TUTORIAL, "Tutorial"),
        (DATA, "Data"), (PAGE, "Page"), (TOOL, "Tool"), (NEWS, "News"),
        (BLOG, "Blog"), (BOARD, "Bulletin Board")
    ]

    TOP_LEVEL = {QUESTION, JOB, FORUM, PAGE, BLOG, DATA, TUTORIAL, TOOL, NEWS, BOARD}

    title = models.CharField(max_length=200, null=False)

    # The user that originally created the post.
    author = models.ForeignKey(User, on_delete=models.SET(get_sentinel_user))

    # The project that this post belongs to.
    project = models.ForeignKey(Project, on_delete=models.SET_NULL, null=True)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(User, related_name='editor', null=True,
                                      on_delete=models.SET(get_sentinel_user))

    # Store users contributing to the thread as "tags" to query later.
    thread_users = models.ManyToManyField(User, related_name="thread_users")

    # Indicates the information value of the post.
    rank = models.FloatField(default=0, blank=True)

    # Indicates whether the post has accepted answer.
    has_accepted = models.BooleanField(default=False, blank=True)

    # Post status: open, closed, deleted.
    status = models.IntegerField(choices=STATUS_CHOICES, default=OPEN)

    # The type of the post: question, answer, comment.
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)

    # Number of upvotes for the post
    vote_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of views for the post.
    view_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of answers that a post has.
    reply_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of comments that a post has.
    comment_count = models.IntegerField(default=0, blank=True)

    # Bookmark count.
    book_count = models.IntegerField(default=0)

    # How many people follow that thread.
    subs_count = models.IntegerField(default=0)

    # The total numbers of votes for a top-level post.
    thread_votecount = models.IntegerField(default=0)

    # The total number of comments + answers for a thread
    thread_score = models.IntegerField(default=0, blank=True, db_index=True)

    # Date related fields.
    creation_date = models.DateTimeField(auto_now_add=True, db_index=True)
    lastedit_date = models.DateTimeField(auto_now_add=True, db_index=True)

    # Stickiness of the post.
    sticky = models.BooleanField(default=False)

    # This will maintain the ancestor/descendant relationship bewteen posts.
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True, on_delete=models.SET_NULL)

    # This will maintain parent/child relationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children', on_delete=models.SET_NULL)

    # This is the HTML that the user enters.
    content = models.TextField(default='')

    # This is the  HTML that gets displayed.
    html = models.TextField(default='')

    # The tag value is the canonical form of the post's tags
    tag_val = models.CharField(max_length=100, default="", blank=True)

    # The tag set is built from the tag string and used only for fast filtering
    tags = TaggableManager()

    # What site does the post belong to.
    site = models.ForeignKey(Site, null=True, on_delete=models.SET_NULL)

    uid = models.CharField(max_length=32, unique=True)

    def parse_tags(self):
        return util.split_tags(self.tag_val)

    def add_tags(self, text):

        text = text.strip()
        if not text:
            return
       # Sanitize the tag value
        self.tag_val = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
       # Clear old tags
        tag_list = [x.lower() for x in self.parse_tags()]
        self.tags.clear()
        self.tags.add(*tag_list)
        self.save()

    @property
    def as_text(self):
        "Returns the body of the post after stripping the HTML tags"
        text = bleach.clean(self.content, tags=[], attributes=[], styles={}, strip=True)
        return text

    def peek(self, length=300):
        "A short peek at the post"
        return self.as_text[:length]

    def get_title(self):
        if self.status == Post.OPEN:
            return self.title
        else:
            return f"({self.get_status_display()}) {self.title}"

    @property
    def is_open(self):
        return self.status == Post.OPEN

    @property
    def is_comment(self):
        return self.type == Post.COMMENT

    @property
    def age_in_days(self):
        delta = timezone.now() - self.creation_date
        return delta.days

    def get_absolute_url(self):
        return reverse("post_view", kwargs=dict(uid=self.root.uid))

    def update_reply_count(self):
        "This can be used to set the answer count."
        if self.type == Post.ANSWER:
            reply_count = Post.objects.filter(parent=self.parent, type=Post.ANSWER, status=Post.OPEN).count()
            Post.objects.filter(pk=self.parent_id).update(reply_count=reply_count)

    @staticmethod
    def update_post_views(post, request, minutes=settings.POST_VIEW_MINUTES):
        "Views are updated per user session"

        # Extract the IP number from the request.
        ip1 = request.META.get('REMOTE_ADDR', '')
        ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
        # 'localhost' is not a valid ip address.
        ip1 = '' if ip1.lower() == 'localhost' else ip1
        ip2 = '' if ip2.lower() == 'localhost' else ip2
        ip = ip1 or ip2 or '0.0.0.0'

        now = util.now()
        since = now - datetime.timedelta(minutes=minutes)

        # One view per time interval from each IP address.
        if not PostView.objects.filter(ip=ip, post=post, date__gt=since).exists():
            PostView.objects.create(ip=ip, post=post, date=now)
            Post.objects.filter(pk=post.pk).update(view_count=F('view_count') + 1)
        return post

    def save(self, *args, **kwargs):

        self.uid = self.uid or util.get_uuid(8)
        self.lastedit_user = self.lastedit_user or self.author

        # Sanitize the post body.
        self.html = self.html or mistune.markdown(self.content)

        # Must add tags with instance method. This is just for safety.
        self.tag_val = util.strip_tags(self.tag_val)

        # Posts other than a question also carry the same tag
        if self.is_toplevel and self.type != Post.QUESTION:
            required_tag = self.get_type_display().lower()

            if self.tag_val and (required_tag not in self.tag_val.split()):
                self.tag_val += "," + required_tag
            else:
                self.tag_val = required_tag

        super(Post, self).save(*args, **kwargs)

        thread = Post.objects.filter(status=self.OPEN, root=self.root)
        reply_count = thread.exclude(pk=self.parent.pk).filter(type=self.ANSWER).count()
        thread_score = thread.exclude(pk=self.root.pk).count()

        Post.objects.filter(pk=self.parent.pk).update(reply_count=reply_count)
        Post.objects.filter(pk=self.root.pk).update(thread_score=thread_score)

    def __str__(self):
        return "%s: %s (pk=%s)" % (self.get_type_display(), self.title, self.pk)

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL

    @property
    def deleted_class(self):
        return "deleted" if self.status == Post.DELETED else ""

    @property
    def accepted_class(self):
        if self.status == Post.DELETED:
            return "deleted"
        return "accepted" if (self.has_accepted and not self.is_toplevel) else ""

    @property
    def belongs_to_project(self):
        return self.project is not None


class Vote(models.Model):
    # Post statuses.

    UP, DOWN, BOOKMARK, ACCEPT, EMPTY = range(5)
    TYPE_CHOICES = [(UP, "Upvote"), (EMPTY, "Empty"),
                    (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.SET(get_user_model))
    post = models.ForeignKey(Post, related_name='votes', on_delete=models.CASCADE)
    type = models.IntegerField(choices=TYPE_CHOICES, default=EMPTY, db_index=True)
    date = models.DateTimeField(auto_now_add=True, db_index=True)

    uid = models.CharField(max_length=32, unique=True)

    def __str__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(limit=16)
        super(Vote, self).save(*args, **kwargs)


class PostView(models.Model):
    """
    Keeps track of post views based on IP address.
    """
    ip = models.GenericIPAddressField(default='', null=True, blank=True)
    post = models.ForeignKey(Post, related_name="post_views", on_delete=models.CASCADE)
    date = models.DateTimeField(auto_now_add=True)


class Subscription(models.Model):
    "Connects a post to a user"

    NO_MESSAGES, LOCAL_MESSAGE, EMAIL_MESSAGE, DIGEST_MESSAGES = range(4)
    MESSAGING_CHOICES = [
        (NO_MESSAGES, "Not following"),
        (LOCAL_MESSAGE, "Follow using Local Messages"),
        (EMAIL_MESSAGE, "Follow using Emails"),
        ]

    class Meta:
        unique_together = (("user", "post"))

    uid = models.CharField(max_length=32, unique=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, related_name="subs",on_delete=models.CASCADE)
    type = models.IntegerField(choices=MESSAGING_CHOICES, default=LOCAL_MESSAGE)
    date = models.DateTimeField()

    objects = SubscriptionManager()

    def __str__(self):
        return "%s to %s" % (self.user.name, self.post.title)

    def save(self, *args, **kwargs):
        # Set the date to current time if missing.
        self.date = self.date or util.now()
        self.uid = self.uid or util.get_uuid(limit=16)
        super(Subscription, self).save(*args, **kwargs)

    @staticmethod
    def get_sub(post, user):

        sub = Subscription.objects.filter(post=post, user=user)
        return None if user.is_anonymous else sub


@receiver(post_save, sender=Post)
def set_post(sender, instance, created, *args, **kwargs ):

    if created:
        # Set the titles
        if instance.parent and not instance.title:
            instance.title = instance.parent.title

        if instance.parent and instance.parent.type in (Post.ANSWER, Post.COMMENT):
            # Only comments may be added to a parent that is answer or comment.
            instance.type = Post.COMMENT

        if instance.type is None:
            # Set post type if it was left empty.
            instance.type = Post.COMMENT if instance.parent else Post.FORUM

        # This runs only once upon object creation.
        instance.title = instance.parent.title if instance.parent else instance.title
        instance.lastedit_user = instance.author
        instance.status = instance.status or Post.PENDING
        instance.creation_date = instance.creation_date or timezone.now()
        instance.lastedit_date = instance.creation_date

        # Set the timestamps on the parent
        if instance.type == Post.ANSWER:
            instance.parent.lastedit_date = instance.lastedit_date
            instance.parent.lastedit_user = instance.lastedit_user
            instance.parent.save()


@receiver(post_save, sender=Post)
def check_thread_users(sender, instance, created, *args, **kwargs):

    # Add the author to the thread_users of the root.

    obj = instance.root if instance.root else instance
    users = list(Post.objects.get_all(pk=obj.pk).values_list("thread_users", flat=True))
    user_pks = users + [instance.author.pk]
    user_pks = set(user_pks)

    # Clear theard users
    obj.thread_users.clear()

    users = User.objects.filter(pk__in=user_pks)
    obj.thread_users.add(*users)


@receiver(post_save, sender=Post)
def check_root(sender, instance, created, *args, **kwargs):
    "We need to ensure that the parent and root are set on object creation."

    if created:

        if not (instance.root or instance.parent):
            # Neither root or parent are set.
            instance.root = instance.parent = instance

        elif instance.parent:
            # When only the parent is set the root must follow the parent root.
            instance.root = instance.parent.root

        elif instance.root:
            # The root should never be set on creation.
            raise Exception('Root may not be set on creation')

        if instance.parent.type in (Post.ANSWER, Post.COMMENT):
            # Answers and comments may only have comments associated with them.
            instance.type = Post.COMMENT

        assert instance.root and instance.parent

        if not instance.is_toplevel:
            # Title is inherited from top level.
            instance.title = "%s: %s" % (instance.get_type_display()[0], instance.root.title[:80])
            instance.project = instance.root.project

            if instance.type == Post.ANSWER:
                Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

        instance.save()

