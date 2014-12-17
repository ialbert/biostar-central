from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf import settings
from django.template import Context, Template, Library

register = Library()

@register.inclusion_tag('widgets/recent-votes.html')
def render_recent_votes(votes):
    return dict(votes=votes)

