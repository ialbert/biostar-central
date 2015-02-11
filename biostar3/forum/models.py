from __future__ import absolute_import, division, print_function, unicode_literals

from django.db import models
from django.contrib.auth.models import UserManager, AbstractBaseUser, UserManager, PermissionsMixin, Group, GroupManager
from django.contrib.sites.models import Site
from django.conf import settings
from django.core.urlresolvers import reverse
from taggit.managers import TaggableManager

class MyTaggableManager(TaggableManager):
    def get_internal_type(self):
        return 'ManyToManyField'

# Default groups.
ADMIN_GROUP_NAME = "Admins"
MODERATOR_GROUP_NAME = "Moderators"

MODERATE_USER_PERMISSION = "moderate_user"
MODERATE_POST_PERMISSION = "moderate_post"
BAN_USER_PERMISSION = "ban_user"


# The main user model.
class User(AbstractBaseUser, PermissionsMixin):
    """Represents a registered user of the website."""

    objects = UserManager()

    # Class level constants.
    USER, MODERATOR, ADMIN, BLOG = range(4)
    TYPE_CHOICES = [(USER, "User"), (MODERATOR, "Moderator"), (ADMIN, "Admin"), (BLOG, "Blog")]

    NEW_USER, TRUSTED, SUSPENDED, BANNED = range(4)
    STATUS_CHOICES = ((NEW_USER, 'New User'), (TRUSTED, 'Trusted'), (SUSPENDED, 'Suspended'), (BANNED, 'Banned'))

    class Meta:
        db_table = "users_user"
        permissions = (
            (MODERATE_USER_PERMISSION, "Can moderate a user"),
            (BAN_USER_PERMISSION, "Can ban a user"),
        )

    # Required by Django.
    USERNAME_FIELD = 'email'

    # Default information on every user.
    email = models.EmailField(verbose_name='Email', db_index=True, max_length=255, unique=True, blank=False)
    name = models.CharField(verbose_name='Name', max_length=255, default="Biostar User", blank=False)

    # Fields used by the Django admin.
    # These are different from the Biostar user types even though Django also calls them admin.
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)
    is_staff = models.BooleanField(default=False)

    # This designates a user types and with that permissions.
    type = models.IntegerField(choices=TYPE_CHOICES, default=USER)

    # This designates a user status on whether they are allowed to log in.
    status = models.IntegerField(choices=STATUS_CHOICES, default=NEW_USER)

    # The number of new messages for the user.
    new_messages = models.IntegerField(default=0)

    # The number of badges for the user.
    badges = models.IntegerField(default=0)

    # Activity score computed over a shorter period.
    score = models.IntegerField(default=0)

    # User's recent activity level.
    activity = models.IntegerField(default=0)

    # Display next to a user name.
    flair = models.CharField(verbose_name='Flair', max_length=15, default="")

    # The site this users belongs to.
    site = models.ForeignKey(Site, null=True)

    @property
    def is_suspended(self):
        return (self.status == self.SUSPENDED) or (self.status == self.BANNED)

    def get_absolute_url(self):
        url = reverse("post_view", kwargs=dict(pk=self.id))
        return url

    def get_short_name(self):
        return self.name

class GroupInfo(models.Model):
    "Extra group information"
    author = models.ForeignKey(settings.AUTH_USER_MODEL)
    group = models.OneToOneField(Group)
    creation_date = models.DateTimeField(auto_now_add=True)

class Tag(models.Model):
    """Represents a tag."""

    class Meta:
        db_table = "posts_tag"

    name = models.TextField(max_length=50, db_index=True)
    count = models.IntegerField(default=0)

class Profile(models.Model):
    """
    Maintains information that does not always need to be retreived whe a user is accessed.
    """

    class Meta:
        db_table = "users_profile"

    # Message type selector.
    TYPE_CHOICES = settings.MESSAGING_MAP.items()

    # Digest choices.
    NO_DIGEST, DAILY_DIGEST, WEEKLY_DIGEST, MONTHLY_DIGEST = range(4)

    DIGEST_CHOICES = [(NO_DIGEST, 'Never'), (DAILY_DIGEST, 'Daily'),
                      (WEEKLY_DIGEST, 'Weekly'), (MONTHLY_DIGEST, 'Monthly')]

    # The user that this profile is for.
    user = models.OneToOneField(settings.AUTH_USER_MODEL)

    # Globally unique id used to identify the user in a private feeds
    uuid = models.CharField(null=False, db_index=True, unique=True, max_length=255)

    # The last visit by the user.
    last_login = models.DateTimeField()

    # The last visit by the user.
    date_joined = models.DateTimeField()

    # User provided location.
    location = models.CharField(default="", max_length=255, blank=True)

    # User provided website.
    website = models.URLField(default="", max_length=255, blank=True)

    # Google scholar ID
    scholar = models.CharField(default="", max_length=255, blank=True)

    # Twitter ID
    twitter_id = models.CharField(default="", max_length=255, blank=True)

    # This field is used to select content for the user.
    my_tags = models.TextField(default="", max_length=255, blank=True)

    # Description provided by the user html.
    info = models.TextField(default="", null=True, blank=True)

    # The default notification preferences.
    message_prefs = models.IntegerField(choices=TYPE_CHOICES, default=settings.MESSAGE_DEFAULT)

    # This stores binary flags on users. Their usage is to
    # allow easy subselection of various subsets of users.
    flag = models.IntegerField(default=0)

    # The tag value is the canonical form of the post's tags
    watched_tags = models.CharField(max_length=250, default="", blank=True)

    # The tag set is built from the watch_tag string and is used to trigger actions
    # when a post that matches this tag is set
    tags = models.ManyToManyField(Tag, blank=True, )


class Post(models.Model):
    """Represents a post."""

    class Meta:
        db_table = "posts_post"
        ordering = ["-lastedit_date"]
        permissions = (
            (MODERATE_POST_PERMISSION, "Can moderate a post"),
        )

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

    # Maintains post tags.
    tags = MyTaggableManager()

    title = models.CharField(max_length=250, null=False)

    # The user that originally created the post.
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # The group that this post belongs to.
    group = models.ForeignKey(Group, null=True, blank=True)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='editor')

    # Indicates the information value of the post.
    rank = models.FloatField(default=0, blank=True)

    # Post status: open, closed, deleted.
    status = models.IntegerField(choices=STATUS_CHOICES, default=OPEN)

    # The type of the post: question, answer, comment.
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)

    # Number of upvotes for the post
    vote_count = models.IntegerField(default=0, blank=True, db_index=True)

    # The number of views for the post.
    view_count = models.IntegerField(default=0, blank=True)

    # The number of replies that a post has.
    reply_count = models.IntegerField(default=0, blank=True)

    # The number of comments that a post has.
    comment_count = models.IntegerField(default=0, blank=True)

    # Bookmark count.
    book_count = models.IntegerField(default=0)

    # Indicates indexing is needed.
    changed = models.BooleanField(default=True)

    # How many people follow that thread.
    subs_count = models.IntegerField(default=0)

    # The total score of the thread (used for top level only)
    thread_score = models.IntegerField(default=0, blank=True, db_index=True)

    # Date related fields.
    creation_date = models.DateTimeField(db_index=True)
    lastedit_date = models.DateTimeField(db_index=True)

    # Stickiness of the post.
    sticky = models.BooleanField(default=False, db_index=True)

    # Indicates whether the post has accepted answer.
    has_accepted = models.BooleanField(default=False, blank=True)

    # This will maintain the ancestor/descendant relationship bewteen posts.
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True)

    # This will maintain parent/child replationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children')

    # This is the HTML that the user enters.
    content = models.TextField(default='')

    # This is the  HTML that gets displayed.
    html = models.TextField(default='')

    # The tag value is the canonical form of the post's tags
    tag_val = models.CharField(max_length=100, default="", blank=True)

    # The tag set is built from the tag string and used only for fast filtering
    tag_set = models.ManyToManyField(Tag, blank=True, )

    # What site does the post belong to.
    site = models.ForeignKey(Site, null=True)

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL

    def get_absolute_url(self):
        url = reverse("post_view", kwargs=dict(pk=self.root_id))
        if self.is_toplevel:
            return url
        else:
            return "%s#%s" % (url, self.id)

    def subtype(self):
        "This is a method of the class the Question type is subdivided."
        if self.type == Post.QUESTION:
            if self.has_accepted:
                return "Accepted"
            elif self.reply_count > 0:
                return "Answered"
            return "Unanswered"
        return self.get_type_display()

class PostView(models.Model):
    """Represents a post vote"""

    class Meta:
        db_table = "posts_postview"

    ip = models.GenericIPAddressField(default='', null=True, blank=True)
    post = models.ForeignKey(Post, related_name="post_views")
    date = models.DateTimeField(auto_now=True)


class FederatedContent(models.Model):
    """
    Represents a searchable text sent over from another site.

    This is a content that can be searched. Used when exchanging content
    across sites. Should be a json object that gets used and deserialized only when indexed.
    """
    obj_id = models.IntegerField(default=0, db_index=True)
    domain = models.TextField(default='', null=False, blank=False)
    content = models.TextField(default='', null=False, blank=False)
    changed = models.BooleanField(default=False, blank=True)
    creation_date = models.DateTimeField(db_index=True, auto_now=True)


class Vote(models.Model):
    class Meta:
        db_table = "posts_vote"

    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Upvote"), (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post, related_name='votes')
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    date = models.DateTimeField(db_index=True, auto_now=True)

    def __unicode__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())
