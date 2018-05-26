import bleach
import logging

from django.utils import timezone
from django.db import models
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.sites.models import Site
from django.dispatch import receiver
from django.db.models.signals import post_save, m2m_changed
from django.db.models import F

from biostar.forum import util

User = get_user_model()


logger = logging.getLogger("engine")

def get_sentinel_user():
    return User.objects.get_or_create(username='deleted').first()


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


class Tag(models.Model):
    name = models.TextField(max_length=50, db_index=True)
    count = models.IntegerField(default=0)

    @staticmethod
    def fixcase(name):
        return name.upper() if len(name) == 1 else name.lower()

    def __str__(self):
        return self.name



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
    author = models.ForeignKey(settings.AUTH_USER_MODEL,
                               on_delete=models.SET(get_sentinel_user))

    # The user that edited the post most recently.
    lastedit_user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='editor',
                                      on_delete=models.SET(get_sentinel_user))

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
    root = models.ForeignKey('self', related_name="descendants", null=True, blank=True, on_delete=models.SET_NULL)

    # This will maintain parent/child replationships between posts.
    parent = models.ForeignKey('self', null=True, blank=True, related_name='children', on_delete=models.SET_NULL)

    # This is the HTML that the user enters.
    content = models.TextField(default='')

    # This is the  HTML that gets displayed.
    html = models.TextField(default='')

    # The tag value is the canonical form of the post's tags
    tag_val = models.CharField(max_length=100, default="", blank=True)

    # The tag set is built from the tag string and used only for fast filtering
    tag_set = models.ManyToManyField(Tag, blank=True)

    # What site does the post belong to.
    site = models.ForeignKey(Site, null=True, on_delete=models.SET_NULL)

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
        delta = timezone.now() - self.creation_date
        return delta.days

    def update_reply_count(self):
        "This can be used to set the answer count."
        if self.type == Post.ANSWER:
            reply_count = Post.objects.filter(parent=self.parent, type=Post.ANSWER, status=Post.OPEN).count()
            Post.objects.filter(pk=self.parent_id).update(reply_count=reply_count)

    def delete(self, using=None, keep_parents=False):
        #Collect tag names.
       tag_names = [t.name for t in self.tag_set.all()]

        #While there is a signal to do this it is much faster this way.
       Tag.objects.filter(name__in=tag_names).update(count=F('count') - 1)

        #Remove tags with zero counts.
       Tag.objects.filter(count=0).delete()
       super(Post, self).delete(using=using, keep_parents=keep_parents)

    def save(self, *args, **kwargs):

        # Sanitize the post body.
        self.html = util.parse_html(self.content)

        # Must add tags with instance method. This is just for safety.
        self.tag_val = util.strip_tags(self.tag_val)

        # Posts other than a question also carry the same tag
        if self.is_toplevel and self.type != Post.QUESTION:
            required_tag = self.get_type_display()
            if required_tag not in self.tag_val:
                self.tag_val += "," + required_tag

        # Recompute post reply count
        self.update_reply_count()

        super(Post, self).save(*args, **kwargs)

    def __str__(self):
        return "%s: %s (pk=%s)" % (self.get_type_display(), self.title, self.pk)

    @property
    def is_toplevel(self):
        return self.type in Post.TOP_LEVEL


class Vote(models.Model):
    # Post statuses.
    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Upvote"), (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.SET(get_user_model))
    post = models.ForeignKey(Post, related_name='votes', on_delete=models.CASCADE)
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    date = models.DateTimeField(db_index=True, auto_now=True)

    def __str__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())


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
            #TODO: will this cause recussion max when instance.root = instance.parent = instance
            instance.parent.save()


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

            if instance.type == Post.ANSWER:
                Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

        instance.save()




@receiver(m2m_changed, sender=Post.tag_set.through)
def update_counts(sender, instance, action, pk_set, *args, **kwargs):
    "Applies tag count updates upon post changes"

    if action == 'post_add':
        Tag.objects.filter(pk__in=pk_set).update(count=F('count') + 1)

    if action == 'post_remove':
        Tag.objects.filter(pk__in=pk_set).update(count=F('count') - 1)

    if action == 'pre_clear':
        instance.tag_set.all().update(count=F('count') - 1)
