from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf import settings
from django.template import Context, Template, Library
import hashlib, urllib, random
from datetime import timedelta, datetime
from django.utils.timezone import utc
from biostar3.forum.models import Post


register = Library()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


POST_TYPE_CSS = dict()
@register.simple_tag
def post_type_css(post):
    if post.type == Post.QUESTION:
        if post.has_accepted:
            return "accepted"
        elif post.reply_count > 0:
            return "answered"

        return "unanswered"

    return post.get_type_display()

@register.inclusion_tag('widgets/recent_votes.html')
def recent_votes(votes):
    return dict(votes=votes)


@register.inclusion_tag('widgets/user_link.html')
def user_link(user):
    return dict(user=user)


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


@register.inclusion_tag('widgets/post_user_box.html')
def post_user_box(post):
    return dict(post=post, author=post.author)


@register.inclusion_tag('widgets/tag_bar.html')
def tag_bar(post, show_update=True):
    return dict(post=post, show_update=show_update)


@register.filter
def bignum(number):
    "Reformats numbers with qualifiers as K"
    try:
        value = float(number) / 1000.0
        if value > 10:
            return "%0.fk" % value
        elif value > 1:
            return "%0.1fk" % value
    except ValueError, exc:
        pass
    return str(number)


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


@register.filter
def time_ago(date):
    delta = now() - date
    if delta < timedelta(minutes=1):
        return 'just now'
    elif delta < timedelta(hours=1):
        unit = pluralize(delta.seconds // 60, "minute")
    elif delta < timedelta(days=1):
        unit = pluralize(delta.seconds // 3600, "hour")
    elif delta < timedelta(days=30):
        unit = pluralize(delta.days, "day")
    elif delta < timedelta(days=90):
        unit = pluralize(int(delta.days / 7), "week")
    elif delta < timedelta(days=730):
        unit = pluralize(int(delta.days / 30), "month")
    else:
        diff = delta.days / 365.0
        unit = '%0.1f years' % diff
    return "%s ago" % unit


@register.simple_tag
def gravatar(user, size=80):
    name = user.name
    if user.is_suspended:
        # Removes spammy images for suspended users
        email = 'suspended@biostars.org'
    else:
        email = user.email.encode('utf8')
    hash = hashlib.md5(email).hexdigest(),

    gravatar_url = "https://secure.gravatar.com/avatar/%s?" % hash
    gravatar_url += urllib.urlencode(dict(s=str(size), d='identicon'))

    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, name)
