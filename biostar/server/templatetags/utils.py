from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
import random

register = template.Library()

@register.simple_tag
def rand_num():
    "The purpose of this is to return a random number"
    return " %f " % random.random()

