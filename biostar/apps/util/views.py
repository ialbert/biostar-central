# Create your views here.
import os
from django import forms
from django.views.generic import  UpdateView, DetailView
from django.contrib.flatpages.models import FlatPage
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.conf import settings

import logging

logger = logging.getLogger(__name__)

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))

