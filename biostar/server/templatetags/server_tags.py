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


@register.inclusion_tag('server_tags/navbar.html', takes_context=True)
def navbar(context, user):
    "Renders top navigation bar"
    return {'user': user, 'CATEGORIES': context['CATEGORIES']}
