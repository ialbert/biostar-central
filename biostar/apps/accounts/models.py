from __future__ import print_function, unicode_literals, absolute_import, division
from django.conf import settings
from django.db import models
from django.contrib.auth.models import User
from allauth.exceptions import ImmediateHttpResponse
from allauth.socialaccount.adapter import DefaultSocialAccountAdapter

# Represents a Biostar user across all sites.
class Profile(models.Model):
    pass
