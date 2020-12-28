
from django.db import models


class Manager(models.Manager):

    def get_queryset(self):
        """
        Route all queries to biostar2 database
        """
        return super().get_queryset().using("biostar2")


class UsersProfile(models.Model):
    website = models.CharField(max_length=255)
    info = models.TextField(blank=True, null=True)
    user_id = models.IntegerField(unique=True)
    uuid = models.CharField(unique=True, max_length=255)
    digest_prefs = models.IntegerField()
    opt_in = models.BooleanField()
    twitter_id = models.CharField(max_length=255)
    my_tags = models.TextField()
    flag = models.IntegerField()
    last_login = models.DateTimeField()
    location = models.CharField(max_length=255)
    watched_tags = models.CharField(max_length=100)
    message_prefs = models.IntegerField()
    scholar = models.CharField(max_length=255)
    id = models.IntegerField(primary_key=True)  # AutoField?
    date_joined = models.DateTimeField()

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'users_profile'


class UsersUser(models.Model):
    NEW_USER, TRUSTED, SUSPENDED, BANNED = range(4)
    STATUS_CHOICES = ((NEW_USER, 'New User'), (TRUSTED, 'Trusted'), (SUSPENDED, 'Suspended'), (BANNED, 'Banned'))

    status = models.IntegerField()
    is_staff = models.BooleanField()
    name = models.CharField(max_length=255)
    site_id = models.IntegerField(blank=True, null=True)
    is_active = models.BooleanField()
    id = models.IntegerField(primary_key=True)  # AutoField?
    score = models.IntegerField()
    last_login = models.DateTimeField()
    new_messages = models.IntegerField()
    activity = models.IntegerField()
    is_admin = models.BooleanField()
    password = models.CharField(max_length=128)
    type = models.IntegerField()
    email = models.CharField(unique=True, max_length=255)
    badges = models.IntegerField()
    flair = models.CharField(max_length=15)

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'users_user'

    @property
    def profile(self):
        return UsersProfile.objects.filter(user_id=self.id).first()


class PostsPost(models.Model):
    site_id = models.IntegerField(blank=True, null=True)
    rank = models.FloatField()
    creation_date = models.DateTimeField()
    reply_count = models.IntegerField()
    tag_val = models.CharField(max_length=100)
    id = models.IntegerField(primary_key=True)
    view_count = models.IntegerField()
    thread_score = models.IntegerField()
    title = models.CharField(max_length=200)
    has_accepted = models.BooleanField()
    vote_count = models.IntegerField()
    content = models.TextField()
    parent_id = models.IntegerField(blank=True, null=True)
    comment_count = models.IntegerField()
    html = models.TextField()
    type = models.IntegerField()
    status = models.IntegerField()
    book_count = models.IntegerField()
    root_id = models.IntegerField(blank=True, null=True)
    lastedit_user_id = models.IntegerField()
    sticky = models.BooleanField()
    lastedit_date = models.DateTimeField()
    changed = models.BooleanField()
    subs_count = models.IntegerField()
    author_id = models.IntegerField()

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'posts_post'

    @property
    def parent(self):
        return self.objects.filter(id=self.parent_id).first()

    @property
    def root(self):
        return self.objects.filter(id=self.root_id).first()

    @property
    def author(self):
        return UsersUser.objects.filter(id=self.author_id).first()

    @property
    def lastedit_user(self):
        return UsersUser.objects.filter(id=self.lastedit_user_id).first()


class PostsSubscription(models.Model):
    id = models.IntegerField(primary_key=True)
    user_id = models.IntegerField()
    post_id = models.IntegerField()
    type = models.IntegerField()
    date = models.DateTimeField()

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'posts_subscription'
        unique_together = (('user_id', 'post_id'),)

    @property
    def user(self):
        return UsersUser.objects.filter(id=self.user_id).first()

    @property
    def post(self):
        return PostsPost.objects.filter(id=self.post_id).first()


class PostsVote(models.Model):
    id = models.IntegerField(primary_key=True)
    author_id = models.IntegerField()
    post_id = models.IntegerField()
    type = models.IntegerField()
    date = models.DateTimeField()

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'posts_vote'

    @property
    def author(self):
        return UsersUser.objects.filter(id=self.author_id).first()

    @property
    def post(self):
        return PostsPost.objects.filter(id=self.post_id).first()


class BadgesBadge(models.Model):
    id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=50)
    type = models.IntegerField()
    unique = models.BooleanField()
    count = models.IntegerField()
    desc = models.CharField(max_length=200)
    icon = models.CharField(max_length=250)

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'badges_badge'


class BadgesAward(models.Model):
    badge_id = models.IntegerField()
    user_id = models.IntegerField()
    date = models.DateTimeField()
    context = models.CharField(max_length=1000)

    objects = Manager()

    class Meta:
        managed = False
        db_table = 'badges_award'

    @property
    def badge(self):
        return BadgesBadge.objects.filter(id=self.badge_id).first()
