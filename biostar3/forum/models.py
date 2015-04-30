from __future__ import absolute_import, division, print_function, unicode_literals
import random, hashlib, uuid, logging
import os, urllib, feedparser
from django.db import models, transaction
from django.contrib.auth.models import UserManager, AbstractBaseUser, PermissionsMixin
from django.contrib.sites.models import Site
from django.conf import settings
from django.core.urlresolvers import reverse
from taggit.managers import TaggableManager
from django.utils import timezone
from datetime import datetime, timedelta
from . import html
from biostar3.utils.compat import *

logger = logging.getLogger('biostar')


def get_group_url(group):
    """
    Find the fully qualified url to a group
    """
    site = Site.objects.get(id=settings.SITE_ID)
    if group.domain == settings.DEFAULT_GROUP_DOMAIN:
        netloc = site.domain
    else:
        netloc = site.domain.split(".")
        if settings.SITE_PREPEND_SUBDOMAIN:
            netloc = [group.domain] + netloc
        else:
            netloc[0] = group.domain
        netloc = ".".join(netloc)

    url = "%s://%s" % (settings.SITE_SCHEME, netloc)
    return url


class MyTaggableManager(TaggableManager):
    def get_internal_type(self):
        return 'ManyToManyField'


def right_now():
    return timezone.now()

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))


def make_uuid(size=8):
    u = uuid.uuid4()
    u = str(u)[:size]
    return u

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

    # Moderator types.
    MODERATOR_TYPES = {MODERATOR, ADMIN}

    NEW_USER, TRUSTED, SUSPENDED, BANNED = range(4)
    STATUS_CHOICES = ((NEW_USER, 'New User'), (TRUSTED, 'Trusted'), (SUSPENDED, 'Suspended'), (BANNED, 'Banned'))

    class Meta:
        db_table = "users_user"

    # Required by Django.
    USERNAME_FIELD = 'email'

    # Default information on every user.
    email = models.EmailField(verbose_name='Email', db_index=True, max_length=255, unique=True, blank=False)
    name = models.CharField(verbose_name='Name', max_length=255, default="", blank=False)

    # Unique handle. Allow empty for now.
    handle = models.CharField(verbose_name='Handle', max_length=25, default='', blank=True )

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
    flair = models.CharField(verbose_name='Flair', max_length=255, default="0,0,0,0", blank=True)

    # The site this users belongs to.
    site = models.ForeignKey(Site, null=True, blank=True)

    # The user portrait.
    portrait = models.FileField(upload_to="users/img/%Y/%m/%d", null=True, blank=True)

    @property
    def is_suspended(self):
        return (self.status == self.SUSPENDED) or (self.status == self.BANNED)

    @property
    def scaled_score(self):
        # Turns out people prefer scores to go by 10.
        return self.score * 10

    def get_absolute_url(self):
        # url = reverse("user_view", kwargs=dict(pk=self.id))
        url = reverse("user_view", kwargs=dict(pk=self.id))
        return url

    def get_short_name(self):
        return self.name

    def post_count(self, types):
        return Post.objects.filter(type__in=types, author=self).count()

    def save(self, *args, **kwargs):
        "Actions that need to be performed on every user save."

        if not self.name:
            # When there is no name try to split from email.
            # Triggered mostly on user creation.
            self.name = self.email.split("@")[0]

        super(User, self).save(*args, **kwargs)

    def __unicode__(self):
        return "User: %s (%s)" % (self.id, self.email)


class UserGroup(models.Model):
    """
    Represents a group
    """
    name = models.CharField(max_length=25, unique=True, db_index=True)
    domain = models.CharField(max_length=50, unique=True, db_index=True, default="www")
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name="owners", null=True)
    info = models.TextField(default="default group info")
    html = models.TextField(default="group info turned into html")

    url = models.CharField(default="full url to the group home", max_length=255)
    public = models.BooleanField(default=True)
    visible = models.BooleanField(default=True)
    user_count = models.IntegerField(default=1)
    post_count = models.IntegerField(default=0)
    last_activity = models.DateTimeField()
    creation_date = models.DateTimeField()

    # The logo for the model.
    logo = models.FileField(upload_to="groups/img", null=True, blank=True)

    # A stylysheet that will be loaded for the group.
    css_file = models.FileField(upload_to="groups/css", null=True, blank=True)

    def save(self, *args, **kwargs):
        "Actions that need to be performed on every user save."

        # A few sanity checks performed.
        self.domain = self.domain.lower()
        self.domain = self.domain.strip()
        self.domain = "".join(self.domain.splitlines())
        self.domain = '-'.join(self.domain.split())

        # The group url is set relative to the current site.
        self.url = get_group_url(self)

        # Set the dates.
        self.creation_date = self.creation_date or right_now()
        self.last_activity = self.last_activity or self.creation_date

        self.html = html.sanitize(self.info, user=self.owner)

        super(UserGroup, self).save(*args, **kwargs)

        if self.owner:
            # This is required since there is a chicken and egg problem with
            # the first usergroup.
            if not GroupSub.objects.filter(user=self.owner, usergroup=self):
                GroupSub.objects.create(user=self.owner, usergroup=self)
                GroupPerm.objects.create(user=self.owner, usergroup=self, role=GroupPerm.ADMIN)

    def get_absolute_url(self):
        return reverse("group_redirect", kwargs=dict(pk=self.id))

    def __unicode__(self):
        return "Usergroup: %s" % self.name

def update_usergroups():
    # Move this to a signal. Not done that way because of migrations would be slow.
    for group in UserGroup.objects.all():
        group.user_count = GroupSub.objects.filter(usergroup=group).count()
        group.post_count = Post.objects.filter(usergroup=group).count()
        group.save()

class GroupPerm(models.Model):
    """
    Represents a special permission for a user on a group.
    """

    class Meta:
        unique_together = (("user", "usergroup"),)

    MODERATE, ADMIN = range(2)
    ROLES = [(MODERATE, "Moderator"), (ADMIN, "Admin")]
    user = models.ForeignKey(settings.AUTH_USER_MODEL, db_index=True)
    usergroup = models.ForeignKey(UserGroup)
    role = models.IntegerField(choices=ROLES, default=MODERATE)


class Profile(models.Model):
    """
    Maintains information that does not always need to be retreived whe a user is accessed.
    """

    class Meta:
        db_table = "users_profile"

    # The user that this profile is for.
    user = models.OneToOneField(settings.AUTH_USER_MODEL)

    # Globally unique id used to identify the user in a private feeds
    uuid = models.CharField(null=False, db_index=True, unique=True, max_length=255)

    # The last visit by the user.
    last_login = models.DateTimeField()

    # The date the user joined user.
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

    # Description provided by the user as markdown.
    info = models.TextField(default="", null=True, blank=True, max_length=5000)

    # The info field turned into html.
    html = models.TextField(default="", null=True, blank=True)

    # The default notification preferences.
    message_prefs = models.IntegerField(choices=settings.MESSAGE_CHOICES, default=settings.MESSAGE_DEFAULT)

    # Subscription to daily and weekly digests.
    digest_prefs = models.IntegerField(choices=settings.DIGEST_CHOICES, default=settings.DEFAULT_DIGEST)

    # This stores binary flags on users. Their usage is to
    # allow easy subselection of various subsets of users.
    flag = models.IntegerField(default=0)

    # The tag value is the canonical form of the post's tags
    watched_tags = models.CharField(max_length=250, default="", blank=True)

    # The text input for the shortcuts.
    shortcuts_text = models.TextField(default="", null=True, blank=True, max_length=1000)

    # Json formatted shortcuts.
    shortcuts_json = models.TextField(default="", null=True, blank=True)

    def save(self, *args, **kwargs):
        self.uuid = self.uuid or make_uuid()
        self.info = self.info.strip()
        self.html = html.sanitize(self.info, user=self.user)
        super(Profile, self).save(*args, **kwargs)


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
    # Some types are not used anymore but are kept for compatibility.
    QUESTION, ANSWER, JOB, FORUM, PAGE, BLOG, COMMENT, DATA, TUTORIAL, BOARD, TOOL, NEWS = range(12)

    TYPE_CHOICES = [
        (QUESTION, "Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (TUTORIAL, "Tutorial"),
        (DATA, "Data"), (PAGE, "Page"), (TOOL, "Tool"), (NEWS, "News"),
        (BLOG, "Blog"), (BOARD, "Bulletin Board")
    ]

    # Posts for which the root=self.
    TOP_LEVEL = {QUESTION, JOB, FORUM, DATA, TUTORIAL, TOOL, NEWS, PAGE}

    # The types of posts that show up in queries.
    QUERY_TYPES = {QUESTION, JOB, FORUM, TUTORIAL, TOOL, NEWS}

    # Maintains post tags.
    tags = MyTaggableManager()

    title = models.CharField(max_length=250, null=False)

    # A unique id for the post.
    uuid = models.CharField(max_length=256, null=True)

    # The user that originally created the post.
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # The group that this post belongs to.
    usergroup = models.ForeignKey(UserGroup, null=True, blank=True)

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='editor', null=True)

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
    last_activity = models.DateTimeField(db_index=True)

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

    # What site does the post belong to.
    site = models.ForeignKey(Site, null=True)

    # The post may have an uploaded file.
    file = models.FileField(upload_to="files/%Y/%m/%d", null=True, blank=True)

    def __unicode__(self):
        display = self.get_type_display()
        return "Post (id={}, type={}, title={})".format(self.id, display, self.title)

    def set_reply_count(self):
        reply_count = Post.objects.filter(parent=self, type=Post.ANSWER).count()
        Post.objects.filter(pk=self.id).update(reply_count=reply_count)

    def get_title(self):
        if self.status == Post.OPEN:
            return self.title
        else:
            return "(%s) %s" % (self.get_status_display(), self.title)

    @property
    def can_accept(self):
        return self.type == Post.ANSWER

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL

    def get_absolute_url(self):

        if self.root:
            url = reverse("post_view", kwargs=dict(pk=self.root_id))
        else:
            # This normally should never happen. It may occur if
            # the object is not updated after creating it. It has self referential Foreign keys.
            # One bad post may make the whole page fail so better protect against it.
            logger.error("missing self.root for id={}".format(self.id))
            url = reverse("post_view", kwargs=dict(pk=self.id))

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

    def update_reply_count(self):
        "This can be used to set the answer count."
        if self.type == Post.ANSWER:
            reply_count = Post.objects.filter(parent=self.parent, type=Post.ANSWER, status=Post.OPEN).count()
            Post.objects.filter(pk=self.parent_id).update(reply_count=reply_count)

    def save(self, *args, **kwargs):
        "Actions that need to be performed on every post save."

        # Attempt to sensibly set the post type if not specified.
        if self.type is None:
            if self.parent.type in (Post.ANSWER, Post.COMMENT):
                self.type = Post.COMMENT
            elif self.parent.type in Post.TOP_LEVEL:
                self.type = Post.ANSWER
            else:
                self.type = Post.QUESTION

        # Set the types of the object
        if self.is_toplevel:
            self.root = self.parent = self
        else:

            if not self.parent:
                raise Exception("non toplevel posts must have a parent")

            if not self.parent.root:
                raise Exception("post parent root not set")

            self.root = self.parent.root

            # Titles are set automatically for content.
            self.title = "%s: %s" % (self.get_type_display()[0], self.root.title)

        # Set the title for the object.
        if not self.title:
            self.title = "%s: %s" % (self.get_type_display()[0], self.parent.title)

        # Update the reply count.
        self.set_reply_count()

        # Remove whitespace from the title.
        self.title = self.title.strip()

        self.creation_date = self.creation_date or right_now()
        self.lastedit_date = self.lastedit_date or self.creation_date
        self.last_activity = self.last_activity or self.lastedit_date
        self.lastedit_user = self.lastedit_user or self.author
        self.html = html.sanitize(self.content, user=self.lastedit_user)

        self.changed = True
        super(Post, self).save(*args, **kwargs)


class PostView(models.Model):
    """Represents a post vote"""

    class Meta:
        db_table = "posts_postview"

    ip = models.GenericIPAddressField(default='', null=True, blank=True)
    post = models.ForeignKey(Post, related_name="post_views")
    date = models.DateTimeField(auto_now=True)


class ReplyToken(models.Model):
    """
    Connects a user and a post to a unique token.
    Sending back the token identifies
    both the user and the post that they are replying to.
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post)
    token = models.CharField(max_length=256)
    date = models.DateTimeField(auto_created=True)

    def save(self, *args, **kwargs):
        self.token = self.token or make_uuid()
        self.date = self.date or right_now()
        super(ReplyToken, self).save(*args, **kwargs)


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


class GroupSub(models.Model):
    """
    Keeps track of the subscription of a user to a group.
    """

    class Meta:
        unique_together = (("user", "usergroup"),)

    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    usergroup = models.ForeignKey(UserGroup)
    type = models.IntegerField(choices=settings.SUBSCRIPTION_CHOICES, default=settings.SUBSCRIPTION_DEFAULT)

    def __unicode__(self):
        return "GroupSub of %s to %s" % (self.user_id, self.usergroup_id)


class PostSub(models.Model):
    """
    Keeps track of subscriptions by users to posts.
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post)
    type = models.IntegerField(choices=settings.SUBSCRIPTION_CHOICES, default=settings.SUBSCRIPTION_DEFAULT)

    class Meta:
        unique_together = (("user", "post"),)

    def __unicode__(self):
        return "PostSub: %s, %s: %s" % (self.user_id, self.post_id, self.get_pref_display())


class Vote(models.Model):
    class Meta:
        db_table = "posts_vote"

    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Upvote"), (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL, db_index=True)
    post = models.ForeignKey(Post, related_name='votes', db_index=True)
    type = models.IntegerField(choices=TYPE_CHOICES)
    unread = models.BooleanField(default=True)
    date = models.DateTimeField(db_index=True)

    def save(self, **kwargs):
        self.date = self.date or right_now()
        super(Vote, self).save(**kwargs)

    def __unicode__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())


class MessageBody(models.Model):
    """
    A message body generated by a user.
    """
    MAX_LEN = 500
    author = models.ForeignKey(settings.AUTH_USER_MODEL, related_name="message_bodies")
    subject = models.CharField(max_length=MAX_LEN)
    html = models.TextField()
    content = models.TextField()
    date = models.DateTimeField()

    def save(self, **kwargs):
        self.subject = self.subject[:self.MAX_LEN]
        self.date = self.date or right_now()
        super(MessageBody, self).save(**kwargs)


class Message(models.Model):
    """
    Connects users to messages
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='messages')
    body = models.ForeignKey(MessageBody, related_name='messages')
    unread = models.BooleanField(default=True)
    date = models.DateTimeField(db_index=True, null=True)

    def save(self, **kwargs):
        self.date = self.date or right_now()
        super(Message, self).save(**kwargs)

    def __unicode__(self):
        return "Message for user %s" % self.user_id


class Blog(models.Model):
    """
    Represents a blog
    """
    title = models.CharField(verbose_name='Blog Name', max_length=255, default="", blank=False)
    desc = models.TextField(default='', blank=True)
    feed = models.URLField()
    link = models.URLField()
    active = models.BooleanField(default=True)
    list_order = models.IntegerField(default=0)
    usergroup = models.ForeignKey(UserGroup, null=True)

    class Meta:
        db_table = "planet_blog"

    @property
    def fname(self):
        fname = abspath(settings.PLANET_DIR, '%s.xml' % self.id)
        return fname

    def parse(self):
        try:
            doc = feedparser.parse(self.fname)
        except Exception as exc:
            logger.error("error %s parsing blog %s", (exc, self.id))
            doc = None
        return doc

    def download(self):
        try:
            text = urlopen(self.feed).read()
            stream = open(self.fname, 'wt')
            stream.write(text)
            stream.close()
        except Exception as exc:
            logger.error("error %s downloading %s", (exc, self.feed))

    def __unicode__(self):
        return self.title

    def get_absolute_url(self):
        return self.link



class BlogPost(models.Model):
    "Represents an entry of a Blog"

    # The blog that generated the entry
    blog = models.ForeignKey(Blog)

    # A unique id for this entry
    uid = models.CharField(max_length=200, default="", null=False)

    # The title of the entry
    title = models.CharField(max_length=200, null=False)

    # The content of the feed
    content = models.TextField(default='', max_length=20000)

    # Sanitized HTML
    html = models.TextField(default='')

    # Date related fields.
    creation_date = models.DateTimeField(db_index=True)

    # Date at which the post has been inserted into the database
    insert_date = models.DateTimeField(db_index=True, null=True)

    # Has the entry been published
    published = models.BooleanField(default=False)

    # The link to the entry
    link = models.URLField()


    class Meta:
        db_table = "planet_blogpost"

    @property
    def get_title(self):
        return u"BLOG: %s" % self.title

    @property
    def is_toplevel(self):
        return True

    def get_absolute_url(self):
        return self.link

    def save(self, *args, **kwargs):
        if not self.id:
            # Set the date to current time if missing.
            self.insert_date = self.insert_date or right_now()
        self.html = html.sanitize(self.content, user=None)
        super(BlogPost, self).save(*args, **kwargs)

    def __unicode__(self):
        return self.title


class Badge(models.Model):
    USER, POST = range(2)
    TARGET_CHOICES = [(USER, "User badge"), (POST, "Post badge")]

    BRONZE, SILVER, GOLD = range(3)
    STYLE_CHOICES = ((BRONZE, 'Bronze'), (SILVER, 'Silver'), (GOLD, 'Gold'))

    # This is used to uniquely identify badges when initializing them.
    uuid = models.CharField(unique=True, default='', blank=True, max_length=100)

    # The name of the badge.
    name = models.CharField(max_length=50)

    # The description of the badge.
    desc = models.CharField(max_length=200, default='')

    # What type of content is it awarded for.
    type = models.IntegerField(choices=TARGET_CHOICES, default=USER)

    # The rarity of the badge.
    style = models.IntegerField(choices=STYLE_CHOICES, default=BRONZE)

    # Unique badges may be earned only once
    unique = models.BooleanField(default=False)

    # Total number of times awarded
    count = models.IntegerField(default=0)

    # The icon to display for the badge.
    icon = models.CharField(default='<i class="fa fa-asterisk"></i>', max_length=250)

    def get_absolute_url(self):
        url = reverse("badge-details", kwargs=dict(pk=self.id))
        return url

    def __unicode__(self):
        return self.name

    def save(self, *args, **kwargs):
        self.uuid = self.uuid or make_uuid()
        super(Badge, self).save(*args, **kwargs)


class Award(models.Model):
    '''
    A badge being awarded to a user.
    '''
    badge = models.ForeignKey(Badge)
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # Keep track of the post if applicable.
    post = models.ForeignKey(Post, null=True, blank=True)
    date = models.DateTimeField()
    context = models.CharField(max_length=1000, default='')

    def save(self, *args, **kwargs):
        User.objects.filter(pk=self.user_id).update(score=models.F('score') + 10)
        super(Award, self).save(*args, **kwargs)

class FlatPage(models.Model):
    # Flatpage is a post that connects to a slug.
    slug = models.SlugField(default='slug', db_index=True)
    post = models.ForeignKey(Post)


def compute_user_score(user, start=None, end=None):
    """
    Computes the score
    """
    now = right_now()

    if start is None:
        start = user.profile.date_joined
    else:
        start = now - timedelta(days=start)

    if end is None:
        end = now
    else:
        end = now - timedelta(days=end)

    vote_count = Vote.objects.filter(post__author=user,
                                     type=Vote.UP, date__gt=start, date__lt=end).count()
    book_count = Vote.objects.filter(post__author=user, type=Vote.BOOKMARK,
                                     date__gt=start, date__lt=end).count()
    post_count = Post.objects.filter(author=user,
                                     creation_date__gt=start, creation_date__lt=end).count()
    answer_count = Post.objects.filter(author=user, type=Post.ANSWER,
                                       creation_date__gt=start, creation_date__lt=end).count()

    score = vote_count * 10 + book_count * 20 + post_count * 1 + answer_count * 5

    return int(score)


def compute_flair(user):
    s1 = compute_user_score(user, start=360, end=90)
    s2 = compute_user_score(user, start=270, end=90)
    s3 = compute_user_score(user, start=180, end=90)
    s4 = compute_user_score(user, start=90, end=0)

    return "{},{},{},{}".format(s1, s2, s3, s4)

