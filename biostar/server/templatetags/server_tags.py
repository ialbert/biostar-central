from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
from biostar.apps.posts.models import Post
import random

register = template.Library()


@register.simple_tag
def rand_num():
    "The purpose of this is to return a random number"
    return " %f " % random.random()


@register.simple_tag
def active(x, y):
    # Create the active class css
    return 'active' if x == y else ''

@register.simple_tag
def boxclass(post):
    # Create the active class css
    if post.type == Post.QUESTION:
        if post.reply_count == 0:
            return "unanswered"
        else:
            return "default"
    else:
        return post.get_type_display()

@register.inclusion_tag('server_tags/navbar.html', takes_context=True)
def navbar(context, user):
    "Renders top navigation bar"
    return {'user': user, 'TOPICS': context['TOPICS'], 'topic': context.get('topic')}


@register.inclusion_tag('server_tags/pagebar.html', takes_context=True)
def pagebar(context):
    "Renders a paging bar"
    return context

@register.inclusion_tag('server_tags/userlink.html')
def userlink(user):
    "User display in a link"
    return {'user': user}


@register.inclusion_tag('server_tags/searchbar.html')
def searchbar():
    "Displays search bar"
    return {}

