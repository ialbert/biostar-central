
from django.db import models


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

    class Meta:
        managed = False
        db_table = 'users_profile'


class UsersUser(models.Model):
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

    class Meta:
        managed = False
        db_table = 'users_user'

    @property
    def profile(self):
        return UsersProfile.objects.using("biostar2").filter(user_id=self.id).first()