from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime, reversion
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.utils.timezone import utc
from django.utils.translation import ugettext_lazy as _
from django.core.urlresolvers import reverse
import bleach
from django.db.models import Q, F
from django.core.exceptions import ObjectDoesNotExist
from biostar import const

# HTML sanitization parameters.
ALLOWED_TAGS = bleach.ALLOWED_TAGS + settings.ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + settings.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(settings.ALLOWED_ATTRIBUTES)

logger = logging.getLogger(__name__)


class Tag(models.Model):
    name = models.TextField(max_length=50, db_index=True)
    count = models.IntegerField(default=0)

    @staticmethod
    def update_counts(sender, instance, action, pk_set, *args, **kwargs):
        "Applies tag count updates upon post changes"

        if action == 'post_add':
            Tag.objects.filter(pk__in=pk_set).update(count=F('count') + 1)

        if action == 'post_delete':
            Tag.objects.filter(pk__in=pk_set).update(count=F('count') - 1)

        if action == 'pre_clear':
            instance.tag_set.all().update(count=F('count') - 1)


class TagAdmin(admin.ModelAdmin):
    list_display = ('name', 'count')
    search_fields = ['name']


admin.site.register(Tag, TagAdmin)


class PostManager(models.Manager):
    def my_bookmarks(self, user):
        query = self.filter(votes__author=user, votes__type=Vote.BOOKMARK).select_related("author")
        query = query.prefetch_related("tag_set")
        return query

    def my_posts(self, user):
        query = self.filter(author=user).select_related("author")
        query = query.prefetch_related("tag_set")
        return query

    def tag_search(self, text):
        "Performs a query by one or more + separated tags"
        text = ''.join(text.split())
        include, exclude = [], []
        for term in text.split('+'):
            if term.endswith("!"):
                exclude.append(term[:-1])
            else:
                include.append(term)

        if include:
            query = self.filter(type__in=Post.TOP_LEVEL, tag_set__name__in=include).exclude(
                tag_set__name__in=exclude).defer('content')
        else:
            query = self.filter(type__in=Post.TOP_LEVEL).exclude(tag_set__name__in=exclude).defer('content')

        query = query.select_related("author").prefetch_related("tag_set").distinct()

        return query

    def get_thread(self, root):
        # Populate the object to build a tree that contains all posts in the thread.
        query = self.filter(root=root).select_related("root author").order_by("type", "-has_accepted", "-vote_count")
        return query

    def top_level(self, user):
        "Returns posts based on a user type"
        if user.is_moderator:
            query = self.filter(type__in=Post.TOP_LEVEL).defer("content")
        else:
            query = self.filter(type__in=Post.TOP_LEVEL, status=Post.OPEN).defer("content")

        return query.select_related("author").prefetch_related("tag_set")


class Post(models.Model):
    "Represents a post in Biostar"

    objects = PostManager()

    #tags = TaggableManager()

    # Post statuses.
    PENDING, OPEN, CLOSED, DELETED = range(4)
    STATUS_CHOICES = [(PENDING, "Pending"), (OPEN, "Open"), (CLOSED, "Closed"), (DELETED, "Deleted")]

    # Question types. Answers should be listed before comments.
    QUESTION, ANSWER, JOB, FORUM, PAGE, BLOG, COMMENT = range(7)
    TYPE_CHOICES = [
        (QUESTION, "Question"), (ANSWER, "Answer"), (COMMENT, "Comment"),
        (JOB, "Job"), (FORUM, "Forum"), (PAGE, "Page"), (BLOG, "Blog"),
    ]

    TOP_LEVEL = set((QUESTION, JOB, FORUM, PAGE, BLOG))

    title = models.CharField(max_length=140, null=False)

    # The user that originally created the post.
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

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

    # This is the sanitized HTML for display.
    content = models.TextField(default='')

    # The tag value is the canonical form of the post's tags
    tag_val = models.CharField(max_length=100, default="", blank=True)

    # The tag set is built from the tag string and used only for fast filtering
    tag_set = models.ManyToManyField(Tag, blank=True, )

    def parse_tags(self):
        names = self.tag_val.split()
        return map(unicode, names)

    def add_tags(self, text):
        text = text.strip()
        # Sanitize the tag value
        self.tag_val = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
        # Clear old tags
        self.tag_set.clear()
        tags = [Tag.objects.get_or_create(name=name)[0] for name in self.parse_tags()]
        self.tag_set.add(*tags)
        self.save()

    def peek(self, length=300):
        text = bleach.clean(self.content, tags=[], attributes=[], styles={}, strip=True)
        text = text[:length]
        return text

    def get_title(self):
        if self.status == Post.OPEN:
            return self.title
        else:
            return "%s [%s]" % (self.title, self.get_status_display())

    @property
    def is_open(self):
        return self.status == Post.OPEN

    @property
    def age_in_days(self):
        delta = const.now() - self.creation_date
        return delta.days


    def save(self, *args, **kwargs):

        # Sanitize the post body.
        self.content = bleach.clean(self.content, tags=ALLOWED_TAGS,
                                    attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)

        # Must add tags with instance method. This is just for safety.
        self.tag_val = bleach.clean(self.tag_val, tags=[], attributes=[], styles={}, strip=True)

        if not self.id:

            # Set the titles
            if self.parent and not self.title:
                self.title = self.parent.title

            if self.parent and self.parent.type in (Post.ANSWER, Post.COMMENT):
                # Only comments may be added to a parent that is answer or comment.
                self.type = Post.COMMENT

            if self.type is None:
                # Set post type if it was left empty.
                self.type = self.COMMENT if self.parent else self.FORUM

            # This runs only once upon object creation.
            self.title = self.parent.title if self.parent else self.title
            self.lastedit_user = self.author
            self.status = self.status or Post.PENDING
            self.creation_date = datetime.datetime.utcnow().replace(tzinfo=utc)
            self.lastedit_date = self.creation_date

        super(Post, self).save(*args, **kwargs)

    def __unicode__(self):
        return "%s: %s (id=%s)" % (self.get_type_display(), self.title, self.id)

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL

    def get_absolute_url(self):
        "A blog will redirect to the original post"
        #if self.url:
        #    return self.url
        url = reverse("post-details", kwargs=dict(pk=self.root_id))
        return url if self.is_toplevel else "%s#%s" % (url, self.id)


    @staticmethod
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

            if not instance.is_toplevel:
                # Title is inherited from top level.
                instance.title = "%s: %s" % (instance.get_type_display()[0], instance.root.title[:80])

            assert instance.root and instance.parent

            if not instance.is_toplevel:
                Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

            instance.save()


# Posts will have revisions.
#reversion.register(Post)

# Revision admin setup.
class PostAdmin(reversion.VersionAdmin):
    list_display = ('title', 'type', 'author')
    fieldsets = (
        (None, {'fields': ('title',)}),
        ('Attributes', {'fields': ('type', 'status', 'sticky',)}),
        ('Content', {'fields': ('tags', 'html', )}),
    )
    search_fields = ('title', 'author__name')


admin.site.register(Post, PostAdmin)


class Vote(models.Model):
    # Post statuses.
    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Upvote"), (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post, related_name='votes')
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    date = models.DateTimeField(db_index=True, auto_now=True)

    def __unicode__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())


class SubscriptionManager(models.Manager):
    def get_subs(self, post):
        "Returns all suscriptions for a post"
        return self.filter(post=post.root).select_related("user")

# This contains the notification types.
from biostar.const import LOCAL_MESSAGE, MESSAGING_TYPE_CHOICES


class Subscription(models.Model):
    "Connects a post to a user"

    class Meta:
        unique_together = (("user", "post"),)

    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name=_("User"), db_index=True)
    post = models.ForeignKey(Post, verbose_name=_("Post"), related_name="subs", db_index=True)
    type = models.IntegerField(choices=MESSAGING_TYPE_CHOICES, default=LOCAL_MESSAGE, db_index=True)
    date = models.DateTimeField(_("Date"), db_index=True)

    objects = SubscriptionManager()

    def __unicode__(self):
        return "%s to %s" % (self.user.name, self.post.title)

    def save(self, *args, **kwargs):

        if not self.id:
            # Set the date to current time if missing.
            self.date = self.date or const.now()

        super(Subscription, self).save(*args, **kwargs)


    @staticmethod
    def get_sub(post, user):

        if user.is_authenticated():
            try:
                return Subscription.objects.get(post=post, user=user)
            except ObjectDoesNotExist, exc:
                return None

        return None

    @staticmethod
    def create(sender, instance, created, *args, **kwargs):
        "Creates a subscription of a user to a post"
        user = instance.author
        root = instance.root
        if Subscription.objects.filter(post=root, user=user).count() == 0:
            sub = Subscription(post=root, user=user, type=user.profile.message_prefs)
            sub.date = datetime.datetime.utcnow().replace(tzinfo=utc)
            sub.save()
            # Increase the subscription count of the root.
            Post.objects.filter(pk=root.id).update(subs_count=F('subs_count') + 1)

    @staticmethod
    def finalize_delete(sender, instance, *args, **kwargs):
        # Decrease the subscription count of the post.
        Post.objects.filter(pk=instance.post.root_id).update(subs_count=F('subs_count') - 1)


# Admin interface for subscriptions
class SubscriptionAdmin(admin.ModelAdmin):
    search_fields = ('user__name', 'user__email')
    list_select_related = ["user", "post"]


admin.site.register(Subscription, SubscriptionAdmin)

# Data signals
from django.db.models.signals import post_save, post_delete, m2m_changed

post_save.connect(Post.check_root, sender=Post)
post_save.connect(Subscription.create, sender=Post, dispatch_uid="create_subs")
post_delete.connect(Subscription.finalize_delete, sender=Subscription, dispatch_uid="delete_subs")
m2m_changed.connect(Tag.update_counts, sender=Post.tag_set.through)

