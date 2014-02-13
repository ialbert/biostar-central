from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
from biostar.apps.posts.models import Post
import random, hashlib, urllib
from datetime import datetime, timedelta
from django.utils.timezone import utc
from django import template

register = template.Library()


@register.simple_tag
def rand_num():
    "The purpose of this is to return a random number"
    return " %f " % random.random()


@register.simple_tag
def gravatar(user, size=80):
    name = user.name
    email = user.email.encode('utf8')
    hash = hashlib.md5(email).hexdigest(),

    gravatar_url = "http://www.gravatar.com/avatar/%s?" % hash
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


@register.simple_tag
def time_ago(post):
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
    return "%s" % unit


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


@register.inclusion_tag('server_tags/nav_bar.html', takes_context=True)
def nav_bar(context, user):
    "Renders top navigation bar"
    return context

@register.inclusion_tag('server_tags/page_bar.html', takes_context=True)
def page_bar(context):
    "Renders a paging bar"
    return context


@register.inclusion_tag('server_tags/post_body.html', takes_context=True)
def post_body(context, post, user, tree):
    "Renders the post body"
    return dict(post=post, user=user, tree=tree, request=context['request'])


@register.inclusion_tag('server_tags/search_bar.html')
def search_bar():
    "Displays search bar"
    return {}


@register.inclusion_tag('server_tags/post_actions.html')
def post_actions(post, user):
    "Renders post actions"
    return dict(post=post, user=user)


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