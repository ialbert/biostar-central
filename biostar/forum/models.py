from django.db import models
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.sites.models import Site

from django.db.models import F

User = get_user_model()




def get_sentinel_user():
    return User.objects.get_or_create(username='deleted').first()



class Vote(models.Model):
    # Post statuses.
    UP, DOWN, BOOKMARK, ACCEPT = range(4)
    TYPE_CHOICES = [(UP, "Upvote"), (DOWN, "DownVote"), (BOOKMARK, "Bookmark"), (ACCEPT, "Accept")]

    author = models.ForeignKey(settings.AUTH_USER_MODEL)
    post = models.ForeignKey(Post, related_name='votes')
    type = models.IntegerField(choices=TYPE_CHOICES, db_index=True)
    date = models.DateTimeField(db_index=True, auto_now=True)

    def __str__(self):
        return u"Vote: %s, %s, %s" % (self.post_id, self.author_id, self.get_type_display())


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

    def __str__(self):
        return self.name


class Post(models.Model):
    "Represents a post in Biostar"

    #objects = PostManager()

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
    site = models.ForeignKey(Site, null=True, on_delete=models.CASCADE)



