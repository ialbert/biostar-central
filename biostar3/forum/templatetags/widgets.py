from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf import settings
from django.template import Context, Template, Library
import hashlib, urllib, random
from datetime import timedelta, datetime
from django.utils.timezone import utc
from biostar3.forum.models import Post
from django.template import loader
from django.core.context_processors import csrf
from django import forms

register = Library()

class SortForm(forms.Form):
    fields = forms.TypedChoiceField(coerce=int, choices=settings.POST_SORT_CHOICES)

class TimeLimit(forms.Form):
    fields = forms.TypedChoiceField(coerce=int, choices=settings.TIME_LIMIT_CHOICES)


def now():
    return datetime.utcnow().replace(tzinfo=utc)

@register.simple_tag
def scoreline(user, size=5):
    random.seed(user.id)
    values = map(lambda x: random.randint(0, 10), range(size))
    #values = [0,1,0,0,0]
    values = map(str, values)
    random.seed()
    return ",".join(values)

@register.simple_tag
def cachebuster():
    value = random.random()
    param = "?x=%f" % value
    return param

@register.simple_tag
def group_logo_img(request):
    group = request.group
    print (group.logo.url)
    return group.logo.url

@register.inclusion_tag('widgets/recent_votes.html')
def recent_votes(votes):
    return dict(votes=votes)

@register.inclusion_tag('widgets/visual_editor.html')
def visual_editor(user, content=''):
    return dict(content=content, user=user)

@register.inclusion_tag('widgets/post_unit.html', takes_context=True)
def post_unit(context, post, comments):
    request = context['request']
    return dict(post=post, comments=comments, request=request)

@register.inclusion_tag('widgets/user_link.html')
def user_link(user):
    return dict(user=user)

@register.filter
def on_value(value):
    "Turn a truth value into an on/off string"
    return "on" if value else 'off'

@register.filter
def nicer_value(value):
    "Show the value if exists or an empty string"
    return value if value else ''

@register.inclusion_tag('widgets/search_bar.html', takes_context=True)
def search_bar(context, page=None, action='search', placeholder="Search"):
    q = context.get('q', '')
    return dict(page=page, q=q, action=action, placeholder=placeholder)

@register.inclusion_tag('widgets/page_bar.html', takes_context=True)
def page_bar(context, page=None):
    return dict(page=page)

@register.inclusion_tag('widgets/action_bar.html')
def action_bar(post, label="ADD COMMENT"):
    return dict(post=post, label=label)

@register.inclusion_tag('widgets/update_bar.html')
def update_bar(post):
    return dict(post=post)

@register.inclusion_tag('widgets/message_bar.html', takes_context=True)
def message_bar(context):
    messages = context.get("messages", '')
    return dict(messages=messages)

@register.inclusion_tag('widgets/tag_bar.html')
def tag_bar(post):
    return dict(post=post)

@register.inclusion_tag('widgets/nav_bar.html', takes_context=True)
def nav_bar(context, user):
    request = context['request']
    return dict(user=user,request=request)

@register.inclusion_tag('widgets/user_bar.html', takes_context=True)
def user_bar(context, user):
    return dict(user=user, context=context)

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


# this contains the body of each comment
COMMENT_TEMPLATE = 'widgets/comment.html'
COMMENT_BODY = loader.get_template(COMMENT_TEMPLATE)
START_TAG = '<div class="indent">'
END_TAG = "</div>"

@register.simple_tag
def render_comments(request, post, tree):
    global COMMENT_BODY, COMMENT_TEMPLATE
    if settings.DEBUG:
        # Reload the template on each request
        COMMENT_BODY = loader.get_template(COMMENT_TEMPLATE)
    if post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree)
    else:
        text = ''
    return text

def traverse_comments(request, post, tree):
    "Traverses the tree and generates the page"
    global COMMENT_BODY, START_TAG, END_TAG

    def traverse(node, collect=[]):
        cont = Context({"post": node, 'user': request.user})
        #cont.update(csrf(request))
        html = COMMENT_BODY.render(cont)
        collect.append(START_TAG)
        collect.append(html)
        for child in tree.get(node.id, []):
            traverse(child, collect=collect)
        collect.append(END_TAG)

    collect = []
    for node in tree.get(post.id, []):
        traverse(node, collect=collect)
    return '\n'.join(collect)