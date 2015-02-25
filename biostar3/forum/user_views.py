from __future__ import absolute_import, division, print_function, unicode_literals

# Python modules.
from collections import OrderedDict, defaultdict

# Django specific modules.
from django.shortcuts import render, redirect, render_to_response
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, RedirectView
from django.template import RequestContext
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from braces.views import LoginRequiredMixin
# Biostar specific local modules.
from . import models, query, search
from . import post_views

# Get custom user model.
User = get_user_model()

class MeView(LoginRequiredMixin, RedirectView):
    """
    A shortcut that redirects user to their account.
    """
    def get_redirect_url(self, *args, **kwargs):
        return self.request.user.get_absolute_url()