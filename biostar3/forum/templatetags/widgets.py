from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf import settings
from django.template import Context, Template, Library

register = Library()

@register.inclusion_tag('widgets/recent_votes.html')
def recent_votes(votes):
    return dict(votes=votes)


@register.inclusion_tag('widgets/page_bar.html', takes_context=True)
def page_bar(context):
    if context.get('is_paginated'):
        page = context['page_obj']
    else:
        page = None
    return dict(page=page, context=context)

@register.inclusion_tag('widgets/search_bar.html', takes_context=True)
def search_bar(context, action='search'):
    q = context.get('q', '')
    posts = context.get('posts', '')
    return dict(q=q, posts=posts, action=action)