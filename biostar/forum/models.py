import logging

import bleach
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.sites.models import Site
from django.db import models
from django.db.models import F
from django.db.models.signals import post_save
from django.dispatch import receiver
from taggit.managers import TaggableManager

from biostar.utils.shortcuts import reverse
from . import util

User = get_user_model()

# The maximum length in characters for a typical name and text field.
MAX_NAME_LEN = 256
MAX_FIELD_LEN = 1024
MAX_TEXT_LEN = 10000
MAX_LOG_LEN = 20 * MAX_TEXT_LEN

logger = logging.getLogger("engine")


class SubscriptionManager(models.Manager):
    def get_subs(self, post):
        "Returns all subscriptions for a post, exclude the "
        return self.filter(post=post.root).select_related("user", "user__profile")


class Post(models.Model):
    "Represents a post in a forum"

    # Post statuses.
    PENDING, OPEN, CLOSED, DELETED = range(4)
    STATUS_CHOICES = [(PENDING, "Pending"), (OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types. Answers should be listed before comments.
    QUESTION, ANSWER, JOB, FORUM, PAGE, BLOG, COMMENT, DATA, TUTORIAL, BOARD, TOOL, NEWS = range(12)

    # Valid post types.
    TYPE_CHOICES = [
        (QUESTION, "Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (TUTORIAL, "Tutorial"),
        (DATA, "Data"), (PAGE, "Page"), (TOOL, "Tool"), (NEWS, "News"),
        (BLOG, "Blog"), (BOARD, "Bulletin Board")
    ]
    TOP_LEVEL = {QUESTION, JOB, FORUM, BLOG, TUTORIAL, TOOL, NEWS}

    # Possile spam states.
    SPAM, NOT_SPAM, DEFAULT = range(3)
    SPAM_CHOICES = [(SPAM, "Spam"), (NOT_SPAM, "Not spam"), (DEFAULT, "Default")]

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

    # The user that last contributed to the thread.
    last_contributor = models.ForeignKey(User, related_name='contributor', null=True,
                                         on_delete=models.CASCADE)

    # Store users contributing to the thread as "tags" to query later.
    thread_users = models.ManyToManyField(User, related_name="thread_users")

    # Indicates the information value of the post.
    rank = models.FloatField(default=0, blank=True, db_index=True)

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

    # Unique id for the post.
    uid = models.CharField(max_length=32, unique=True, db_index=True)

    # Spam labeling.
    spam = models.IntegerField(choices=SPAM_CHOICES, default=DEFAULT)

    def parse_tags(self):
        return util.split_tags(self.tag_val)

    def add_tags(self, text):

        text = text.strip()
        if not text:
            return
        # Sanitize the tag value
        # Clear old tags
        # self.tags.clear()
        # self.tags.add(*tag_list)

    @property
    def as_text(self):
        "Returns the body of the post after stripping the HTML tags"
        text = bleach.clean(self.content, tags=[], attributes={}, styles=[], strip=True)
        return text

    @property
    def get_votecount(self):

        if self.is_toplevel:
            return self.thread_votecount
        return self.vote_count

    @property
    def get_reply_count(self):
        return self.reply_count

    @property
    def is_open(self):
        return self.status == Post.OPEN

    @property
    def is_deleted(self):
        return self.status == Post.DELETED

    @property
    def has_accepted(self):
        return bool(self.accept_count)

    @property
    def is_comment(self):
        return self.type == Post.COMMENT

    def get_absolute_url(self):
        return reverse("post_view", kwargs=dict(uid=self.root.uid))

    def save(self, *args, **kwargs):

        # Needs to be imported here to avoid circular imports.
        from biostar.utils import markdown

        self.lastedit_user = self.lastedit_user or self.author
        self.creation_date = self.creation_date or util.now()
        self.lastedit_date = self.lastedit_date or self.creation_date
        self.last_contributor = self.lastedit_user

        # Sanitize the post body.
        self.html = markdown.parse(self.content)

        # Set the rank
        self.rank = self.lastedit_date.timestamp()

        # Must add tags with instance method. This is just for safety.
        self.tag_val = util.strip_tags(self.tag_val)

        self.creation_date = self.creation_date or util.now()

        self.lastedit_date = self.lastedit_date or self.creation_date

        if self.type == Post.ANSWER:
            Post.objects.filter(uid=self.parent.uid).update(lastedit_date=self.lastedit_date,
                                                            lastedit_user=self.lastedit_user)

        # This will trigger the signals
        super(Post, self).save(*args, **kwargs)

    def __str__(self):
        return "%s: %s (pk=%s)" % (self.get_type_display(), self.title, self.pk)

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL

    @property
    def css(self):
        # Used to simplify CSS rendering.
        return self.get_status_display()

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

    # NO_MESSAGES, LOCAL_MESSAGE, EMAIL_MESSAGE, DIGEST_MESSAGES = range(4)
    LOCAL_MESSAGE, EMAIL_MESSAGE, NO_MESSAGES, DEFAULT_MESSAGES, ALL_MESSAGES = range(5)

    MESSAGING_CHOICES = [
        (DEFAULT_MESSAGES, "Default to Local Messages"),
        (LOCAL_MESSAGE, "Local messages"),
        (EMAIL_MESSAGE, "Email for every new post added to current one."),
        (ALL_MESSAGES, "Email for every new thread (mailing list mode)")
    ]

    class Meta:
        unique_together = (("user", "post"))

    uid = models.CharField(max_length=32, unique=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, related_name="subs", on_delete=models.CASCADE)
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
        self.uid = self.uid or util.get_uuid(limit=4)
        super(Badge, self).save(*args, **kwargs)


class Award(models.Model):
    '''
    A badge being awarded to a user.Cannot be ManyToManyField
    because some may be earned multiple times
    '''
    badge = models.ForeignKey(Badge, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    post = models.ForeignKey(Post, null=True, on_delete=models.SET_NULL)
    date = models.DateTimeField()
    # context = models.CharField(max_length=1000, default='')
    uid = models.CharField(max_length=32, unique=True)

    def save(self, *args, **kwargs):
        # Set the date to current time if missing.
        self.uid = self.uid or util.get_uuid(limit=16)
        super(Award, self).save(*args, **kwargs)


@receiver(post_save, sender=Post)
def complete_post(sender, instance, created, *args, **kwargs):
    # Determine the root of the post.
    root = instance.root if instance.root is not None else instance

    # Make current user first in the list of contributors.
    root.thread_users.remove(instance.author)
    root.thread_users.add(instance.author)

    if created:
        # Make the Uid user friendly
        instance.uid = instance.uid or f"p{instance.pk}"

        # Set the titles
        if instance.parent and not instance.title:
            instance.title = instance.parent.title

        # Only comments may be added to a parent that is answer or comment.
        if instance.parent and instance.parent.type in (Post.ANSWER, Post.COMMENT):
            instance.type = Post.COMMENT

        # Set post type if it was left empty.
        if instance.type is None:
            instance.type = Post.COMMENT if instance.parent else Post.FORUM

        # This runs only once upon object creation.
        instance.title = instance.parent.title if instance.parent else instance.title
        instance.lastedit_user = instance.author
        instance.last_contributor = instance.author
        instance.status = instance.status or Post.PENDING

        # Default tags
        instance.tag_val = instance.tag_val or "tag1,tag2"

        if instance.parent:
            # When the parent is set the root must follow the parent root.
            instance.root = instance.parent.root
        else:
            # When there is no parent, root and parent are set to itself.
            instance.root = instance.parent = instance

        # Answers and comments may only have comments associated with them.
        if instance.parent.type in (Post.ANSWER, Post.COMMENT):
            instance.type = Post.COMMENT

        # Sanity check.
        assert instance.root and instance.parent

        # Title is inherited from top level.
        if not instance.is_toplevel:
            instance.title = "%s: %s" % (instance.get_type_display()[0], instance.root.title[:80])
        # Update the root answer count
        if instance.type == Post.ANSWER:
            Post.objects.filter(id=instance.root.id).update(answer_count=F("answer_count") + 1)
        # Update root comment count
        if instance.type == Post.COMMENT:
            Post.objects.filter(id=instance.root.id).update(comment_count=F("comment_count") + 1)
        # Update the reply counts.
        thread = Post.objects.filter(status=instance.OPEN, root=instance.root)
        reply_count = thread.exclude(uid=instance.parent.uid).count()
        # Update the instance reply count
        instance.reply_count = reply_count

        # Update the root reply count
        Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

        # Update last contributor to the thread.
        instance.root.last_contributor = instance.last_contributor

        instance.save()
