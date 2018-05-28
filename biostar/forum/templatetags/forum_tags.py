import logging
from datetime import timedelta, datetime
from django.utils.timezone import utc
from django import template


logger = logging.getLogger("engine")

register = template.Library()

def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.inclusion_tag('wideg/post_body.html', takes_context=True)
def post_body(context, post, user, tree):
    "Renders the post body"
    return dict(post=post, user=user, tree=tree, request=context['request'])


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)

@register.filter
def time_ago(date):
    # Rare bug. TODO: Need to investigate why this can happen.
    if not date:
        return ''
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



@register.filter
def bignum(number):
    "Reformats numbers with qualifiers as K"
    try:
        value = float(number) / 1000.0
        if value > 10:
            return "%0.fk" % value
        elif value > 1:
            return "%0.1fk" % value
    except ValueError as exc:
        pass
    return str(number)


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
