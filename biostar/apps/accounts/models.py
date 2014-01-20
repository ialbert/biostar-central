from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django import forms
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.forms import ReadOnlyPasswordHashField
from django.contrib.auth.models import BaseUserManager, AbstractBaseUser, UserManager
from django.utils.timezone import utc

logger = logging.getLogger(__name__)


class User(AbstractBaseUser):
    # Class level constants.
    NEW, MEMBER, MODERATOR, ADMIN = range(4)
    USER_TYPE_CHOICES = [(NEW, "New"), (MEMBER, "Member"), (MODERATOR, "Moderator"), (ADMIN, "Admin")]

    ACTIVE, SUSPENDED, BANNED = range(3)
    USER_STATUS_CHOICES = ( (ACTIVE, 'Active'), (SUSPENDED, 'Suspended'), (BANNED, 'Banned' ))

    # Required by Django.
    USERNAME_FIELD = 'email'

    objects = UserManager()

    # Default information on every user.
    email = models.EmailField(verbose_name='email address', db_index=True, max_length=255, unique=True)
    name = models.CharField(verbose_name='display name', max_length=255, default="", blank=False)

    # Fields used by the Django admin.
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)
    is_staff = models.BooleanField(default=False)

    # This designates a user types and with that permissions.
    type = models.IntegerField(choices=USER_TYPE_CHOICES, default=NEW)

    # This designates a user statuses on whether they are allowed to log in.
    status = models.IntegerField(choices=USER_STATUS_CHOICES, default=ACTIVE)

    # User total reputation.
    reputation = models.IntegerField(default=0)

    # The number of new messages for the user.
    new_messages = models.IntegerField(default=0)

    # The number of badges for the user.
    badges = models.IntegerField(default=0)

    # Activity score computed over a shorter period.
    score = models.IntegerField(default=0)

    def get_full_name(self):
        # The user is identified by their email address
        return self.name or self.email

    def get_short_name(self):
        # The user is identified by their email address
        return self.name or self.email

    def has_perm(self, perm, obj=None):
        "Does the user have a specific permission?"
        # Simplest possible answer: Yes, always
        return True

    def has_module_perms(self, app_label):
        "Does the user have permissions to view the app `app_label`?"
        # Simplest possible answer: Yes, always
        return True

    def save(self, *args, **kwargs):
        "Contains the actions that need to be peformed on every user save."

        if not self.name:
            # Name must always be set.
            self.name = self.email.split("@")[0]
            logger.info("setting name to %s" % self.name)

        super(User, self).save(*args, **kwargs)

    def __unicode__(self):
        return unicode(self.email)


class Profile(models.Model):
    """
    Maintains information that does not always need to be retreived whe a user is accessed.
    """
    user = models.OneToOneField(User)

    # Globally unique id used to identify the user in a private feeds
    uuid = models.TextField(null=False, db_index=True, unique=True)

    # The last visit by the user.
    last_login = models.DateTimeField()

    # The last visit by the user.
    date_joined = models.DateTimeField()

    # User provided location.
    location = models.TextField(default="", null=True, blank=True)

    # User provided website.
    website = models.URLField(default="", null=True, max_length=250, blank=True)

    # This field is used to select content for the user.
    my_tags = models.TextField(default="", null=True, max_length=250, blank=True)

    # Google scholar ID
    scholar = models.TextField(null=True, default='', max_length=50, blank=True)

    # Description provided by the user as markup
    info = models.TextField(default="", null=True, blank=True)

    # The markup rendered as html.
    info_html = models.TextField(default="", null=True, blank=True)

    def save(self, *args, **kwargs):
        # Generate html from the markdown.
        self.info_html = self.info

        if not self.id:
            # This runs only on object creation.
            self.date_joined = datetime.datetime.utcnow().replace(tzinfo=utc)
            self.last_login = self.date_joined

        super(Profile, self).save(*args, **kwargs)


class UserCreationForm(forms.ModelForm):
    """A form for creating new users."""
    password1 = forms.CharField(label='Password', widget=forms.PasswordInput)
    password2 = forms.CharField(label='Password confirmation', widget=forms.PasswordInput)

    class Meta:
        model = User
        fields = ('email', 'name')

    def clean_password2(self):
        # Check that the two password entries match
        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError("Passwords don't match")
        return password2

    def save(self, commit=True):
        # Save the provided password in hashed format
        user = super(UserCreationForm, self).save(commit=False)
        user.set_password(self.cleaned_data["password1"])
        if commit:
            user.save()
        return user

class UserChangeForm(forms.ModelForm):
    """A form for updating users."""
    password = ReadOnlyPasswordHashField()

    class Meta:
        model = User
        fields = ['email', 'password', 'name', 'type', 'is_active', 'is_admin', 'is_staff']

    def clean_password(self):
        # Regardless of what the user provides, return the initial value.
        # This is done here, rather than on the field, because the
        # field does not have access to the initial value
        return self.initial["password"]


class BiostarUserAdmin(UserAdmin):
    # The forms to add and change user instances
    form = UserChangeForm
    add_form = UserCreationForm

    # The fields to be used in displaying the User model.
    # These override the definitions on the base UserAdmin
    # that reference specific fields on auth.User.
    list_display = ('email', 'name', 'type', 'is_admin', 'is_staff')
    list_filter = ('is_admin',)
    fieldsets = (
        (None, {'fields': ('email', 'password')}),
        ('Personal info', {'fields': ('name', 'type')}),
        ('Permissions', {'fields': ('is_admin', 'is_staff')}),
    )
    # add_fieldsets is not a standard ModelAdmin attribute. UserAdmin
    # overrides get_fieldsets to use this attribute when creating a user.
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('email', 'name', 'type', 'password1', 'password2')}
        ),
    )
    search_fields = ('email',)
    ordering = ('email',)
    filter_horizontal = ()

# Register the class in the admin interface.
admin.site.register(User, BiostarUserAdmin)

# Data signals
from django.db.models.signals import post_save

def create_profile(sender, instance, created, *args, **kwargs):
    if created:
        prof = Profile(user=instance)
        prof.save()

post_save.connect(create_profile, sender=User)