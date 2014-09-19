from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
from biostar.apps.posts.models import Post, Tag
from biostar.apps.messages.models import Message
import random, hashlib, urllib
from datetime import datetime, timedelta
from django.utils.timezone import utc
from django import template
from django.core.urlresolvers import reverse
from biostar import const
from biostar.server.views import LATEST

register = template.Library()

@register.simple_tag
def get_count(counts, word):
    num = counts.get(word.lower()) or ''
    return num

@register.simple_tag
def current(request, *urls):
    if request.path in (reverse(url) for url in urls):
        return "active"
    return ''

# The purpose of this is to return a random number
# that makes resources look different and therefore reload
if settings.DEBUG:
    @register.simple_tag
    def rand_num():
        return " %f " % random.random()
else:
    # Turns it off when not in debug mode.
    @register.simple_tag
    def rand_num():
        return "1"

@register.filter
def show_nonzero(value):
    "The purpose of this is to return value or empty"
    return value if value else ''

@register.filter
def bignum(number):
    "Reformats numbers with qualifiers as K"
    try:
        value = float(number)/1000.0
        if value > 10:
            return "%0.fk" % value
        elif value > 1:
            return "%0.1fk" % value
    except ValueError, exc:
        pass
    return str(number)

@register.filter
def on(value):
    "The purpose of this is to return value or empty"
    return "on" if value else 'off'

@register.filter
def latest(value):
    "Attempts to hide parts of the email"
    print "-" * 10, value
    return value if value else "Latest"

@register.filter
def hide_email(value):
    "Attempts to hide parts of the email"
    try:
        addr, host = value.split('@')
        hide = '*' * (len(addr) - 1)
        email = addr[0] + hide + '@' + host
        return email
    except Exception, exc:
        return value

@register.simple_tag
def messages_read(user):
    Message.objects.filter(user=user, unread=True).update(unread=False)
    return ''

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
    gravatar_url += urllib.urlencode({
        's': str(size),
        'd': 'identicon',
    }
    )
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, name)


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


@register.filter
def time_ago(date):
    delta = const.now() - date
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
def last_action(post):
    action = "written"
    return "%s" % action

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

@register.inclusion_tag('server_tags/sidebar_posts.html')
def sidebar_posts(posts):
    return dict(posts=posts)

@register.inclusion_tag('server_tags/sidebar_votes.html')
def sidebar_votes(votes):
    return dict(votes=votes)

@register.inclusion_tag('server_tags/sidebar_users.html')
def sidebar_users(users):
    return dict(users=users)

@register.inclusion_tag('server_tags/sidebar_locations.html')
def sidebar_locations(users):
    return dict(users=users)

@register.inclusion_tag('server_tags/sidebar_awards.html')
def sidebar_awards(awards):
    return dict(awards=awards)

@register.inclusion_tag('server_tags/nav_bar.html', takes_context=True)
def nav_bar(context, user):
    "Renders top navigation bar"
    return context

@register.inclusion_tag('server_tags/page_bar.html', takes_context=True)
def page_bar(context):
    "Renders a paging bar"
    return context

@register.inclusion_tag('server_tags/post_user_box.html')
def post_user_box(user, date):
    "Renders a user box"
    return dict(user=user, date=date)

@register.inclusion_tag('server_tags/user_box.html')
def user_box(user, lastlogin):
    "Renders a user box"
    return dict(user=user, lastlogin=lastlogin)

@register.inclusion_tag('server_tags/page_bar_sort_posts.html', takes_context=True)
def page_bar_sort_posts(context):
    context['sort_fields'] = const.POST_SORT_FIELDS
    context['date_fields'] = const.POST_LIMIT_FIELDS
    "Renders a paging bar"
    return context

@register.inclusion_tag('server_tags/page_bar_sort_users.html', takes_context=True)
def page_bar_sort_users(context):
    context['sort_fields'] = const.USER_SORT_FIELDS
    context['date_fields'] = const.POST_LIMIT_FIELDS
    "Renders a paging bar"
    return context


@register.inclusion_tag('server_tags/post_body.html', takes_context=True)
def post_body(context, post, user, tree):
    "Renders the post body"
    return dict(post=post, user=user, tree=tree, request=context['request'])


@register.inclusion_tag('server_tags/search_bar.html', takes_context=True)
def search_bar(context):
    "Displays search bar"
    return context

@register.inclusion_tag('server_tags/post_count_box.html')
def post_count_box(post, context='', topic=''):
    "Displays the count box for a post row"
    topic = Tag.fixcase(topic)
    topic = topic.split('+')
    if LATEST in topic:
        topic.remove(LATEST)
    return dict(post=post, context=context, topic=topic)

@register.inclusion_tag('server_tags/post_actions.html')
def post_actions(post, user, label="COMMENT"):
    "Renders post actions"
    return dict(post=post, user=user, label=label)


@register.inclusion_tag('server_tags/user_link.html')
def userlink(user):
    "Renders the flair"
    marker = "&bull;"
    if user.is_admin:
        marker = '&diams;&diams;'
    elif user.is_moderator:
        marker = '&diams;'
    return {'user': user, 'marker': marker}

# this contains the body of each comment
COMMENT_TEMPLATE = 'server_tags/comment_body.html'
COMMENT_BODY = template.loader.get_template(COMMENT_TEMPLATE)


@register.simple_tag
def render_comments(request, post, tree):
    global COMMENT_BODY, COMMENT_TEMPLATE
    if settings.DEBUG:
        # reload the template to get changes
        COMMENT_BODY = template.loader.get_template(COMMENT_TEMPLATE)
    if post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree)
    else:
        text = ''
    return text


def traverse_comments(request, post, tree):
    "Traverses the tree and generates the page"
    global COMMENT_BODY

    def traverse(node):
        data = ['<div class="indent">']
        cont = Context({"post": node, 'user': request.user, 'request': request})
        cont.update(csrf(request))
        html = COMMENT_BODY.render(cont)
        data.append(html)
        for child in tree.get(node.id, []):
            data.append(traverse(child))
        data.append("</div>")
        return '\n'.join(data)

    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:
        coll.append(traverse(node))
    return '\n'.join(coll)