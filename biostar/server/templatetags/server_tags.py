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
    return {'user': user, 'TOPICS': context['TOPICS']}

@register.inclusion_tag('server_tags/pagebar.html')
def pagebar(objs):
    "Renders a paging bar"
    return {'objs': objs}

@register.inclusion_tag('server_tags/userlink.html')
def userlink(user):
    "User display in a link"
    return {'user': user}

@register.inclusion_tag('server_tags/searchbar.html')
def searchbar():
    "Displays search bar"
    return {}