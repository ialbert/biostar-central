from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
from biostar.apps.posts.models import Post
import random, hashlib, urllib
from datetime import datetime, timedelta
from django.utils.timezone import utc

register = template.Library()


@register.simple_tag
def rand_num():
    "The purpose of this is to return a random number"
    return " %f " % random.random()


@register.simple_tag
def gravatar(user, size=80):
    username = user.profile.display_name
    useremail = user.email.encode('utf8')
    hash = hashlib.md5(useremail).hexdigest(),

    gravatar_url = "http://www.gravatar.com/avatar/%s?" % hash
    gravatar_url += urllib.urlencode({
        's': str(size),
        'd': 'identicon',
    }
    )
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, username)


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


@register.simple_tag
def action_time_ago(post):
    now = datetime.utcnow().replace(tzinfo=utc)
    delta = now - post.lastedit_date
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

    action = post.get_update_type_display().lower()
    return "%s %s" % (action, unit)


@register.simple_tag
def active(x, y):
    # Create the active class css
    x, y = x or '', y or ''
    return 'active' if x.lower() == y.lower() else ''


@register.simple_tag
def boxclass(post):
    # Create the css class for each row
    if post.has_accepted:
        style = "accepted"
    elif post.reply_count > 0:
        style = "answered"
    elif post.comment_count > 0:
        style = "commented"
    else:
        style = "unanswered"
    return style


@register.inclusion_tag('server_tags/navbar.html', takes_context=True)
def navbar(context, user):
    "Renders top navigation bar"
    return {'user': user, 'TOPICS': context['TOPICS'], 'topic': context.get('topic')}


@register.inclusion_tag('server_tags/pagebar.html', takes_context=True)
def pagebar(context):
    "Renders a paging bar"
    return context


@register.inclusion_tag('server_tags/searchbar.html')
def searchbar():
    "Displays search bar"
    return {}


@register.inclusion_tag('server_tags/userlink.html')
def userlink(user):
    "Renders the flair"
    marker = "&bull;"
    if user.is_admin:
        marker = '&diams;&diams;'
    elif user.is_moderator:
        marker = '&diams;'
    return {'user': user, 'marker': marker}
