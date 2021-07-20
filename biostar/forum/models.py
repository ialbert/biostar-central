import logging
from datetime import timedelta
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.sites.models import Site
from django.core.cache import cache
from django.core.cache.utils import make_template_fragment_key
from django.db import models
from django.db.models import F
from django.db.models import Q
from django.shortcuts import reverse
from taggit.managers import TaggableManager
from urllib.parse import urlparse

from biostar.utils import helpers
from biostar.accounts.models import Profile
from biostar.planet.models import BlogPost
from . import util

User = get_user_model()

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN

logger = logging.getLogger("engine")


class PostManager(models.Manager):

    def valid_posts(self, u=None, **kwargs):
        """
        Returns posts that are not closed or marked as spam.
        """
        query = super().get_queryset().filter(**kwargs)

        # Exclude posts with no root or parents
        query = query.exclude(Q(root=None) | Q(parent=None))

        # Moderators get to see all posts by default.
        if u and u.is_authenticated and u.profile.is_moderator:
            return query

        # Filter for open posts.
        query = query.filter(status=Post.OPEN, root__status=Post.OPEN)

        # Mark every spam post as 'closed'
        # This doubles the query time.
        query = query.exclude(models.Q(spam=Post.SPAM) | models.Q(root__spam=Post.SPAM))

        return query

    def old(self, **kwargs):
        """
        Return posts that were transferred over from an older verion of biostars
        """
        query = super().get_queryset().exclude(uid__contains='p').filter(**kwargs)
        return query


def delete_fragment_cache(key, *params):
    """
    Drops a template fragment cache.
    """
    key = make_template_fragment_key(key, params)
    cache.delete(key)


def delete_post_cache(post):
    """
    Drops both post specific template fragment caches.
    """
    delete_fragment_cache("post", True, post.uid)
    delete_fragment_cache("post", False, post.uid)
    if post.root:
        delete_fragment_cache("post", True, post.root.uid)
        delete_fragment_cache("post", False, post.root.uid)


class Post(models.Model):
    "Represents a post in a forum"

    # Post statuses.
    PENDING, OPEN, OFFTOPIC, CLOSED, DELETED = range(5)
    STATUS_CHOICES = [(PENDING, "Pending"), (OPEN, "Open"), (OFFTOPIC, "Off topic"), (CLOSED, "Closed"),
                      (DELETED, "Deleted")]

    # Question types. Answers should be listed before comments.
    QUESTION, ANSWER, JOB, FORUM, PAGE, BLOG, COMMENT, DATA, TUTORIAL, BOARD, TOOL, NEWS, HERALD = range(13)

    # Valid post types.
    TYPE_CHOICES = [
        (QUESTION, "Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (TUTORIAL, "Tutorial"),
        (DATA, "Data"), (PAGE, "Page"), (TOOL, "Tool"), (NEWS, "News"),
        (BLOG, "Blog"), (BOARD, "Bulletin Board"), (HERALD, "Herald")
    ]
    TOP_LEVEL = {QUESTION, JOB, FORUM, BLOG, TUTORIAL, TOOL, NEWS, HERALD}

    # Possible spam states.
    SPAM, NOT_SPAM, DEFAULT = range(3)

    # SPAM, NOT_SPAM, DEFAULT, SUSPECT = range(4)
    SPAM_CHOICES = [(SPAM, "Spam"), (NOT_SPAM, "Not spam"), (DEFAULT, "Default")]
    # Spam labeling.
    spam = models.IntegerField(choices=SPAM_CHOICES, default=DEFAULT, db_index=True)

    # Spam score stores relative likely hood this post is spam.
    spam_score = models.FloatField(default=0)

    # Post status: open, closed, deleted.
    status = models.IntegerField(choices=STATUS_CHOICES, default=OPEN, db_index=True)

    # The type of the post: question, answer, comment.
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)

    # Post title.
    title = models.CharField(max_length=200, null=False, db_index=True)

    # The user that originally created the post.
    author = models.ForeignKey(User, on_delete=models.CASCADE)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(User, related_name='editor', null=True,
                                      on_delete=models.CASCADE)

    # Store users contributing to the thread as "tags" to more_like_this later.
    thread_users = models.ManyToManyField(User, related_name="thread_users")

    # Indicates the information value of the post.
    rank = models.FloatField(default=0, blank=True, db_index=True)

    # This post has been indexed by the search engine.
    indexed = models.BooleanField(default=False)

    # Show that post is top level
    is_toplevel = models.BooleanField(default=False, db_index=True)

    # Indicates whether the post has accepted answer.
    answer_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of accepted answers.
    accept_count = models.IntegerField(default=0, blank=True)

    # The number of replies for  thread.
    reply_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of comments that a post has.
    comment_count = models.IntegerField(default=0, blank=True)

    # Number of upvotes for the post
    vote_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The total numbers of votes for a top-level post.
    thread_votecount = models.IntegerField(default=0, db_index=True)

    # The number of views for the post.
    view_count = models.IntegerField(default=0, blank=True, db_index=True)

    # Bookmark count.
    book_count = models.IntegerField(default=0)

    # How many people follow that thread.
    subs_count = models.IntegerField(default=0)

    # Post creation date.
    creation_date = models.DateTimeField(db_index=True)

    # Post last edit date.
    lastedit_date = models.DateTimeField(db_index=True)

    # Sticky posts go on top.
    sticky = models.BooleanField(default=False)

    # This will maintain the ancestor/descendant relationship bewteen posts.
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True, on_delete=models.CASCADE)

    # This will maintain parent/child relationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children', on_delete=models.CASCADE)

    # This is the text that the user enters.
    content = models.TextField(default='')

    # This is the  HTML that gets displayed.
    html = models.TextField(default='')

    # The tag value is the canonical form of the post's tags
    tag_val = models.CharField(max_length=100, default="", blank=True)

    # The tag set is built from the tag string and used only for fast filtering
    tags = TaggableManager()

    # What site does the post belong to.
    site = models.ForeignKey(Site, null=True,  on_delete=models.CASCADE)

    # Unique id for the post.
    uid = models.CharField(max_length=32, unique=True, db_index=True)

    # This post has been indexed by the search engine.
    has_diff = models.BooleanField(default=False)

    objects = PostManager()

    def parse_tags(self):
        return [tag.lower() for tag in self.tag_val.split(",") if tag]

    @property
    def get_votecount(self):

        if self.is_toplevel:
            return self.thread_votecount
        return self.vote_count

    def title_prefix(self):

        prefix = ""
        if self.is_spam:
            prefix = "Spam:"
        elif not self.is_open or not self.is_question:
            prefix = f"{self.get_type_display()}:" if self.is_open else f"{self.get_status_display()}:"

        return prefix

    def is_herald(self):
        return self.type == self.HERALD

    @property
    def is_open(self):
        return self.status == Post.OPEN and not self.is_spam

    def recompute_scores(self):
        # Recompute answers count
        if self.type == Post.ANSWER:
            answer_count = Post.objects.valid_posts(root=self.root, type=Post.ANSWER).count()
            Post.objects.filter(pk=self.parent_id).update(answer_count=answer_count)

        reply_count = Post.objects.valid_posts(root=self.root).exclude(pk=self.root.pk).count()

        Post.objects.filter(pk=self.root.id).update(reply_count=reply_count)

    def json_data(self):
        data = {
            'id': self.id,
            'uid': self.uid,
            'title': self.title,
            'type': self.get_type_display(),
            'type_id': self.type,
            'creation_date': util.datetime_to_iso(self.creation_date),
            'lastedit_date': util.datetime_to_iso(self.lastedit_date),
            'author_uid': self.author.profile.uid,
            'lastedit_user_uid': self.lastedit_user.profile.uid,
            'author': self.author.profile.name,
            'status': self.get_status_display(),
            'status_id': self.status,
            'thread_score': self.thread_votecount,
            'rank': self.rank,
            'vote_count': self.vote_count,
            'view_count': self.view_count,
            'reply_count': self.reply_count,
            'comment_count': self.comment_count,
            'book_count': self.book_count,
            'subs_count': self.subs_count,
            'answer_count': self.root.reply_count,
            'has_accepted': self.has_accepted,
            'parent_id': self.parent.id,
            'root_id': self.root_id,
            'xhtml': self.html,
            'content': self.content,
            'tag_val': self.tag_val,
            'url': f'{settings.PROTOCOL}://{settings.SITE_DOMAIN}{self.get_absolute_url()}',
        }
        return data

    @property
    def is_question(self):
        return self.type == Post.QUESTION

    @property
    def is_job(self):
        return self.type == Post.JOB

    @property
    def is_deleted(self):
        return self.status == Post.DELETED

    @property
    def not_spam(self):
        return self.spam == Post.NOT_SPAM

    @property
    def has_accepted(self):
        return bool(self.accept_count)

    def num_lines(self, offset=0):
        """
        Return number of lines in post content
        """
        return len(self.content.split("\n")) + offset

    @property
    def is_spam(self):
        return self.spam == self.SPAM

    @property
    def is_comment(self):
        return self.type == Post.COMMENT

    @property
    def is_answer(self):
        return self.type == Post.ANSWER

    def get_absolute_url(self):
        url = reverse("post_view", kwargs=dict(uid=self.root.uid))
        return url if self.is_toplevel else "%s#%s" % (url, self.uid)

    def high_spam_score(self, threshold=None):
        threshold = threshold or settings.SPAM_THRESHOLD
        return (self.spam_score > threshold) or self.is_spam or self.author.profile.low_rep

    def save(self, *args, **kwargs):

        # Needs to be imported here to avoid circular imports.
        from biostar.forum import markdown

        self.lastedit_user = self.lastedit_user or self.author

        self.creation_date = self.creation_date or util.now()
        self.lastedit_date = self.lastedit_date or util.now()

        # Sanitize the post body.
        self.html = markdown.parse(self.content, post=self, clean=True, escape=False)
        self.tag_val = self.tag_val.replace(' ', '')
        # Default tags
        self.tag_val = self.tag_val or "tag1,tag2"
        # Set the top level state of the post.
        self.is_toplevel = self.type in Post.TOP_LEVEL

        # Drop the cached fragment
        delete_post_cache(self)

        # This will trigger the signals
        super(Post, self).save(*args, **kwargs)

    def __str__(self):
        return "%s: %s (pk=%s)" % (self.get_type_display(), self.title, self.pk)

    def update_parent_counts(self):
        """
        Update the counts for the parent and root
        """

        descendants = Post.objects.filter(root=self.root).exclude(Q(pk=self.root.pk) | Q(status=Post.DELETED)
                                                                  | Q(spam=Post.SPAM))
        answer_count = descendants.filter(type=Post.ANSWER).count()
        comment_count = descendants.filter(type=Post.COMMENT).count()
        reply_count = descendants.count()
        # Update the root reply, answer, and comment counts.
        Post.objects.filter(pk=self.root.pk).update(reply_count=reply_count, answer_count=answer_count,
                                                    comment_count=comment_count)

        children = Post.objects.filter(parent=self.parent).exclude(pk=self.parent.pk)
        com_count = children.filter(type=Post.COMMENT).count()

        # Update parent reply, answer, and comment counts.
        Post.objects.filter(pk=self.parent.pk, is_toplevel=False).update(comment_count=com_count, answer_count=0,
                                                                         reply_count=children.count())

    @property
    def css(self):
        # Used to simplify CSS rendering.
        status = self.get_status_display()
        if self.is_spam:
            return "spam"

        return f"{status}".lower()

    @property
    def accepted_class(self):
        if self.status == Post.DELETED:
            return "deleted"
        if self.has_accepted and not self.is_toplevel:
            return "accepted"

        return ""

    @property
    def age_in_days(self):
        delta = util.now() - self.creation_date
        return delta.days


class Vote(models.Model):
    # Post statuses.

    UP, DOWN, BOOKMARK, ACCEPT, EMPTY = range(5)

    TYPE_CHOICES = [(UP, "Upvote"), (EMPTY, "Empty"),
                    (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, related_name='votes', on_delete=models.CASCADE)
    type = models.IntegerField(choices=TYPE_CHOICES, default=EMPTY, db_index=True)
    date = models.DateTimeField(db_index=True)

    def __str__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())

    def save(self, *args, **kwargs):
        # self.uid = self.uid or f"v{util.get_uuid(limit=5)}"
        self.date = self.date or util.now()
        super(Vote, self).save(*args, **kwargs)

    @property
    def uid(self):
        return


class PostView(models.Model):
    """
    Keeps track of post views based on IP address.
    """
    ip = models.GenericIPAddressField(default='', null=True, blank=True)
    post = models.ForeignKey(Post, related_name="post_views", on_delete=models.CASCADE)
    date = models.DateTimeField(auto_now_add=True)


def update_post_views(post, request, timeout=settings.POST_VIEW_TIMEOUT):
    """
    Views are updated per interval.
    """

    # Get the ip.
    ip = helpers.get_ip(request)

    # Keys go by IP and post ip.
    cache_key = f"{ip}-{post.id}"

    # Found hit no need to increment the views
    if cache.get(cache_key):
        return

    # Insert a new view into database.
    PostView.objects.create(ip=ip, post=post)

    # Separately increment post view.
    Post.objects.filter(id=post.id).update(view_count=F('view_count') + 1)

    # Set the cache.
    cache.set(cache_key, 1, timeout)

    # Drop the post related cache for logged in users.
    if request.user.is_authenticated:
        delete_post_cache(post)

    return post


class Subscription(models.Model):
    "Connects a post to a user"

    LOCAL_MESSAGE, EMAIL_MESSAGE, NO_MESSAGES = range(3)
    SUB_CHOICES = [(LOCAL_MESSAGE, "Local messages"), (EMAIL_MESSAGE, "Email message"), (NO_MESSAGES, "Not subscribed")]
    TYPE_MAP = {Profile.NO_MESSAGES: NO_MESSAGES,
                Profile.EMAIL_MESSAGE: EMAIL_MESSAGE,
                Profile.LOCAL_MESSAGE: LOCAL_MESSAGE,
                Profile.DEFAULT_MESSAGES: LOCAL_MESSAGE}

    class Meta:
        unique_together = (("user", "post"))

    user = models.ForeignKey(User, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, related_name="subs", on_delete=models.CASCADE)
    type = models.IntegerField(choices=SUB_CHOICES, null=True, default=LOCAL_MESSAGE)
    date = models.DateTimeField()

    def __str__(self):
        return f"{self.user.profile.name} to {self.post.title}"

    def save(self, *args, **kwargs):
        # Set the date to current time if missing.
        self.date = self.date or util.now()
        # self.uid = self.uid or util.get_uuid(limit=16)

        if self.type is None:
            self.type = self.TYPE_MAP.get(self.user.profile.message_prefs, self.NO_MESSAGES)

        super(Subscription, self).save(*args, **kwargs)

    def profile_type_mapper(self):
        type_map = {Profile.NO_MESSAGES: self.NO_MESSAGES,
                    Profile.EMAIL_MESSAGE: self.EMAIL_MESSAGE,
                    Profile.LOCAL_MESSAGE: self.LOCAL_MESSAGE,
                    Profile.DEFAULT_MESSAGES: self.LOCAL_MESSAGE}
        return type_map

    @staticmethod
    def get_sub(post, user):
        sub = Subscription.objects.filter(post=post, user=user).first()
        return None if user.is_anonymous else sub

    @property
    def uid(self):
        return self.pk


class SharedLink(models.Model):

    # User submitting the herald
    author = models.ForeignKey(User, on_delete=models.CASCADE)

    # Sets title
    title = models.CharField(max_length=200, blank=True, default='')

    # User that accepted/decline this submission.
    editor = models.ForeignKey(User, related_name='herald_editor', on_delete=models.SET_NULL, null=True)

    # URL of the given herald_list
    url = models.URLField(max_length=MAX_TEXT_LEN)

    # Text ( markdown ) description and html representation.
    text = models.TextField(max_length=MAX_TEXT_LEN, blank=True, default='')

    # Date this herald_list was created.
    creation_date = models.DateTimeField()
    lastedit_date = models.DateTimeField()

    # Gains a post once published, assumed none until then.
    post = models.ForeignKey(Post, on_delete=models.SET_NULL, null=True)

    SUBMITTED, DECLINED, ACCEPTED, PUBLISHED = range(4)
    CHOICES = [(SUBMITTED, 'Submitted'), (DECLINED, 'Declined'), (PUBLISHED, 'Published'), (ACCEPTED, 'Accepted')]
    status = models.IntegerField(choices=CHOICES, default=SUBMITTED, db_index=True)

    @property
    def domain(self):
        """
        Returns the domain of the url
        """
        try:
            domain = urlparse(self.url).netloc
        except Exception as exc:
            domain = ''

        return domain

    def __str__(self):
        return self.url

    def save(self, *args, **kwargs):
        # Needs to be imported here to avoid circular imports.

        self.creation_date = self.creation_date or util.now()
        self.lastedit_date = self.lastedit_date or self.creation_date or util.now()
        self.title = self.title or ''
        super(SharedLink, self).save(*args, **kwargs)
        return

    @property
    def declined(self):
        return self.status == self.DECLINED

    @property
    def published(self):
        return self.status == self.PUBLISHED

    @property
    def submitted(self):
        return self.status == self.SUBMITTED

    @property
    def accepted(self):
        return self.status == self.ACCEPTED

    @property
    def icon(self):

        if self.accepted:
            return 'green check'
        elif self.published:
            return 'purple book'
        elif self.declined:
            return 'orange times'
        else:
            return 'blue paper plane'


class Badge(models.Model):
    BRONZE, SILVER, GOLD = range(3)
    CHOICES = ((BRONZE, 'Bronze'), (SILVER, 'Silver'), (GOLD, 'Gold'))

    # The name of the badge.
    name = models.CharField(max_length=50)

    # The description of the badge.
    desc = models.CharField(max_length=200, default='')

    # The rarity of the badge.
    type = models.IntegerField(choices=CHOICES, default=BRONZE)

    # The icon to display for the badge.
    icon = models.CharField(default='', max_length=250)

    uid = models.CharField(max_length=32, unique=True)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        # Set the date to current time if missing.
        self.uid = self.uid or util.get_uuid(limit=8)
        super(Badge, self).save(*args, **kwargs)


class Award(models.Model):
    '''
    A badge being awarded to a user.Cannot be ManyToManyField
    because some may be earned multiple times
    '''
    badge = models.ForeignKey(Badge, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, null=True, on_delete=models.CASCADE)
    date = models.DateTimeField()

    # context = models.CharField(max_length=1000, default='')

    def save(self, *args, **kwargs):
        # Set the date to current time if missing.
        self.date = self.date or util.now()
        super(Award, self).save(*args, **kwargs)

    @property
    def uid(self):
        return self.pk


class Diff(models.Model):

    # Initial content state
    diff = models.TextField(default='')

    # Date this change was made.
    created = models.DateTimeField(auto_now_add=True)

    # Post this diff belongs to
    post = models.ForeignKey(Post, on_delete=models.CASCADE)

    # Person who created the diff
    author = models.ForeignKey(User, on_delete=models.CASCADE)

    def save(self, *args, **kwargs):

        self.created = self.created or util.now()

        super(Diff, self).save(*args, **kwargs)

    @property
    def breakline(self):
        diff = self.diff
        diff = diff.replace('\n', '<br>')
        return diff


class Log(models.Model):
    """
    Represents moderation actions
    """
    MODERATE, CREATE, EDIT, LOGIN, LOGOUT, CLASSIFY, DEFAULT = range(7)

    ACTIONS_CHOICES = [
        (MODERATE, "Moderate"),
        (CREATE, "Create"),
        (EDIT, "Edit"),
        (LOGIN, "Login"),
        (LOGOUT, "Logout"),
        (CLASSIFY, "Classify"),
        (DEFAULT, "Default")
    ]

    # User that performed the action.
    user = models.ForeignKey(User, null=True, blank=True, on_delete=models.CASCADE)

    # A potential target user (it may be null)
    target = models.ForeignKey(User, related_name="target", null=True, blank=True, on_delete=models.CASCADE)

    # Post related information goes here (it may be null).
    post = models.ForeignKey(Post, null=True, blank=True, on_delete=models.CASCADE)

    # The IP address associated with the log.
    ipaddr = models.GenericIPAddressField(null=True, blank=True)

    # Actions that the user took.
    action = models.IntegerField(choices=ACTIONS_CHOICES, default=DEFAULT, db_index=True)

    # The logging information.
    text = models.TextField(null=True, blank=True)

    # Date this log was created.
    date = models.DateTimeField()

    def save(self, *args, **kwargs):
        self.date = self.date or util.now()
        super(Log, self).save(*args, **kwargs)

