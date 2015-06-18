from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime, string
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.contrib.sites.models import Site

from django.utils.timezone import utc
from django.utils.translation import ugettext_lazy as _
from django.core.urlresolvers import reverse
import bleach
from django.db.models import Q, F
from django.core.exceptions import ObjectDoesNotExist
from biostar import const
from biostar.apps.util import html
from biostar.apps import util
# HTML sanitization parameters.

logger = logging.getLogger(__name__)

def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)

class Tag(models.Model):
    name = models.TextField(max_length=50, db_index=True)
    count = models.IntegerField(default=0)

    @staticmethod
    def fixcase(name):
        return name.upper() if len(name) == 1 else name.lower()

    @staticmethod
    def update_counts(sender, instance, action, pk_set, *args, **kwargs):
        "Applies tag count updates upon post changes"

        if action == 'post_add':
            Tag.objects.filter(pk__in=pk_set).update(count=F('count') + 1)

        if action == 'post_remove':
            Tag.objects.filter(pk__in=pk_set).update(count=F('count') - 1)

        if action == 'pre_clear':
            instance.tag_set.all().update(count=F('count') - 1)

    def __unicode__(self):
        return self.name

class TagAdmin(admin.ModelAdmin):
    list_display = ('name', 'count')
    search_fields = ['name']


admin.site.register(Tag, TagAdmin)

class PostManager(models.Manager):

    def my_bookmarks(self, user):
        query = self.filter(votes__author=user, votes__type=Vote.BOOKMARK)
        query = query.select_related("root", "author", "lastedit_user")
        query = query.prefetch_related("tag_set")
        return query

    def my_posts(self, target, user):

        # Show all posts for moderators or targets
        if user.is_moderator or user == target:
            query = self.filter(author=target)
        else:
            query = self.filter(author=target).exclude(status=Post.DELETED)

        query = query.select_related("root", "author", "lastedit_user")
        query = query.prefetch_related("tag_set")
        query = query.order_by("-creation_date")
        return query

    def fixcase(self, text):
        return text.upper() if len(text) == 1 else text.lower()

    def tag_search(self, text):
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
            query = self.filter(type__in=Post.TOP_LEVEL, tag_set__name__in=include).exclude(
                tag_set__name__in=exclude)
        else:
            query = self.filter(type__in=Post.TOP_LEVEL).exclude(tag_set__name__in=exclude)

        query = query.filter(status=Post.OPEN)

        # Remove fields that are not used.
        query = query.defer('content', 'html')

        # Get the tags.
        query = query.select_related("root", "author", "lastedit_user").prefetch_related("tag_set").distinct()

        return query

    def get_thread(self, root, user):
        # Populate the object to build a tree that contains all posts in the thread.
        is_moderator = user.is_authenticated() and user.is_moderator
        if is_moderator:
            query = self.filter(root=root).select_related("root", "author", "lastedit_user").order_by("type", "-has_accepted", "-vote_count", "creation_date")
        else:
            query = self.filter(root=root).exclude(status=Post.DELETED).select_related("root", "author", "lastedit_user").order_by("type", "-has_accepted", "-vote_count", "creation_date")

        return query

    def top_level(self, user):
        "Returns posts based on a user type"
        is_moderator = user.is_authenticated() and user.is_moderator
        if is_moderator:
            query = self.filter(type__in=Post.TOP_LEVEL)
        else:
            query = self.filter(type__in=Post.TOP_LEVEL).exclude(status=Post.DELETED)

        return query.select_related("root", "author", "lastedit_user").prefetch_related("tag_set").defer("content", "html")


class Post(models.Model):
    "Represents a post in Biostar"

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

    TOP_LEVEL = set((QUESTION, JOB, FORUM, PAGE, BLOG, DATA, TUTORIAL, TOOL, NEWS, BOARD))

    title = models.CharField(max_length=200, null=False)

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

    def parse_tags(self):
        return util.split_tags(self.tag_val)

    def add_tags(self, text):
        text = text.strip()
        if not text:
            return
        # Sanitize the tag value
        self.tag_val = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
        # Clear old tags
        self.tag_set.clear()
        tags = [Tag.objects.get_or_create(name=name)[0] for name in self.parse_tags()]
        self.tag_set.add(*tags)
        #self.save()

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
            return "(%s) %s" % ( self.get_status_display(), self.title)

    @property
    def is_open(self):
        return self.status == Post.OPEN

    @property
    def age_in_days(self):
        delta = const.now() - self.creation_date
        return delta.days

    def update_reply_count(self):
        "This can be used to set the answer count."
        if self.type == Post.ANSWER:
            reply_count = Post.objects.filter(parent=self.parent, type=Post.ANSWER, status=Post.OPEN).count()
            Post.objects.filter(pk=self.parent_id).update(reply_count=reply_count)

    def delete(self, using=None):
        # Collect tag names.
        tag_names = [t.name for t in self.tag_set.all()]

        # While there is a signal to do this it is much faster this way.
        Tag.objects.filter(name__in=tag_names).update(count=F('count') - 1)

        # Remove tags with zero counts.
        Tag.objects.filter(count=0).delete()
        super(Post, self).delete(using=using)

    def save(self, *args, **kwargs):

        # Sanitize the post body.
        self.html = html.parse_html(self.content)

        # Must add tags with instance method. This is just for safety.
        self.tag_val = html.strip_tags(self.tag_val)

        # Posts other than a question also carry the same tag
        if self.is_toplevel and self.type != Post.QUESTION:
            required_tag = self.get_type_display()
            if required_tag not in self.tag_val:
                self.tag_val += "," + required_tag

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
            self.creation_date = self.creation_date or now()
            self.lastedit_date = self.creation_date

            # Set the timestamps on the parent
            if self.type == Post.ANSWER:
                self.parent.lastedit_date = self.lastedit_date
                self.parent.lastedit_user = self.lastedit_user
                self.parent.save()

        # Recompute post reply count
        self.update_reply_count()

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
    def update_post_views(post, request, minutes=settings.POST_VIEW_MINUTES):
        "Views are updated per user session"

        # Extract the IP number from the request.
        ip1 = request.META.get('REMOTE_ADDR', '')
        ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
        # 'localhost' is not a valid ip address.
        ip1 = '' if ip1.lower() == 'localhost' else ip1
        ip2 = '' if ip2.lower() == 'localhost' else ip2
        ip = ip1 or ip2 or '0.0.0.0'

        now = const.now()
        since = now - datetime.timedelta(minutes=minutes)

        # One view per time interval from each IP address.
        if not PostView.objects.filter(ip=ip, post=post, date__gt=since):
            PostView.objects.create(ip=ip, post=post, date=now)
            Post.objects.filter(id=post.id).update(view_count=F('view_count') + 1)
        return post

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

            assert instance.root and instance.parent

            if not instance.is_toplevel:
                # Title is inherited from top level.
                instance.title = "%s: %s" % (instance.get_type_display()[0], instance.root.title[:80])

                if instance.type == Post.ANSWER:
                    Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

            instance.save()


class ReplyToken(models.Model):
    """
    Connects a user and a post to a unique token. Sending back the token identifies
    both the user and the post that they are replying to.
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post)
    token = models.CharField(max_length=256)
    date = models.DateTimeField(auto_created=True)

    def save(self, *args, **kwargs):
        if not self.id:
            self.token = util.make_uuid()
        super(ReplyToken, self).save(*args, **kwargs)

class ReplyTokenAdmin(admin.ModelAdmin):
    list_display = ('user', 'post', 'token', 'date')
    ordering = ['-date']
    search_fields = ('post__title', 'user__name')

admin.site.register(ReplyToken, ReplyTokenAdmin)


class EmailSub(models.Model):
    """
    Represents an email subscription to the digest digest.
    """
    SUBSCRIBED, UNSUBSCRIBED = 0, 1
    TYPE_CHOICES = [
        (SUBSCRIBED, "Subscribed"), (UNSUBSCRIBED, "Unsubscribed"),

    ]
    email = models.EmailField()
    status = models.IntegerField(choices=TYPE_CHOICES)


class EmailEntry(models.Model):
    """
    Represents an digest digest email entry.
    """
    DRAFT, PENDING, PUBLISHED = 0, 1, 2

    # The email entry may be posted as an entry.
    post = models.ForeignKey(Post, null=True)

    # This is a simplified text content of the Post body.
    text = models.TextField(default='')

    # The data the entry was created at.
    creation_date = models.DateTimeField(auto_now_add=True)

    # The date the email was sent
    sent_at = models.DateTimeField(null=True, blank=True)

    # The date the email was sent
    status = models.IntegerField(choices=((DRAFT, "Draft"), (PUBLISHED, "Published")))


class PostAdmin(admin.ModelAdmin):
    list_display = ('title', 'type', 'author')
    fieldsets = (
        (None, {'fields': ('title',)}),
        ('Attributes', {'fields': ('type', 'status', 'sticky',)}),
        ('Content', {'fields': ('content', )}),
    )
    search_fields = ('title', 'author__name')

admin.site.register(Post, PostAdmin)


class PostView(models.Model):
    """
    Keeps track of post views based on IP address.
    """
    ip = models.GenericIPAddressField(default='', null=True, blank=True)
    post = models.ForeignKey(Post, related_name="post_views")
    date = models.DateTimeField(auto_now=True)


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

class VoteAdmin(admin.ModelAdmin):
    list_display = ('author', 'post', 'type', 'date')
    ordering = ['-date']
    search_fields = ('post__title', 'author__name')


admin.site.register(Vote, VoteAdmin)

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
            sub_type = user.profile.message_prefs
            if sub_type == const.DEFAULT_MESSAGES:
                sub_type = const.EMAIL_MESSAGE if instance.is_toplevel else const.LOCAL_MESSAGE
            sub = Subscription(post=root, user=user, type=sub_type)
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

