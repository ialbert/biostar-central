from __future__ import absolute_import, division, print_function, unicode_literals
import hashlib
import random
from datetime import timedelta, datetime
from django import template

from django.conf import settings
from django.template import Context, Library
from django.utils.timezone import utc
from django.template import loader
from django import forms

from biostar3.forum.models import Vote
from biostar3.utils.compat import *
from biostar3.forum import html

register = Library()


class SortForm(forms.Form):
    fields = forms.TypedChoiceField(coerce=int, choices=settings.POST_SORT_CHOICES)


class TimeLimit(forms.Form):
    fields = forms.TypedChoiceField(coerce=int, choices=settings.TIME_LIMIT_CHOICES)


def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.simple_tag
def scoreline(user):
    return user.flair


@register.simple_tag
def cachebuster():
    value = random.random()
    param = "?x=%f" % value
    return param


@register.simple_tag
def group_logo_img(request):
    group = request.group
    return group.logo.url


@register.inclusion_tag('recent_votes.html')
def recent_votes(votes):
    return dict(votes=votes)


@register.inclusion_tag('recent_users.html')
def recent_users(users):
    return dict(users=users)

@register.inclusion_tag('post_visual_editor.html')
def visual_editor(user, content='', show_upload=False):
    return dict(content=content, user=user, show_upload=show_upload)


@register.inclusion_tag('post_unit.html', takes_context=True)
def post_unit(context, post, comments):
    request = context['request']
    return dict(post=post, comments=comments, request=request)


@register.inclusion_tag('user_link.html')
def user_link(user):
    return dict(user=user)


@register.filter
def on_value(value):
    "Turn a truth value into an on/off string"
    return "on" if value else 'off'


@register.filter
def hide_email(value):
    "Hides parts of an email"
    try:
        addr, host = value.split('@')
        hide = '*' * (len(addr) - 1)
        email = addr[0] + hide + '@' + host
        return email
    except Exception as exc:
        return value


@register.filter
def nicer_value(value):
    "Show the value if exists or an empty string"
    return value if value else ''


@register.inclusion_tag('search_bar.html', takes_context=True)
def search_bar(context, page=None, action='search', placeholder="Search"):
    q = context.get('q', '')
    return dict(page=page, q=q, action=action, placeholder=placeholder, context=context)


@register.inclusion_tag('page_bar.html', takes_context=True)
def page_bar(context, page=None):
    return dict(page=page)


@register.inclusion_tag('post_action_bar.html')
def action_bar(post, label="ADD COMMENT"):
    return dict(post=post, label=label)


@register.inclusion_tag('post_update_bar.html')
def update_bar(post):
    return dict(post=post)


@register.inclusion_tag('message_bar.html', takes_context=True)
def message_bar(context):
    messages = context.get("messages", '')
    return dict(messages=messages)


@register.inclusion_tag('tag_bar.html')
def tag_bar(post):
    return dict(post=post)


@register.inclusion_tag('nav_bar.html', takes_context=True)
def nav_bar(context, user):
    group = context.get('group')
    return dict(user=user, context=context, group=group)


@register.inclusion_tag('user_bar.html', takes_context=True)
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
    except ValueError as exc:
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
def glow(value):
    return "glow" if value else ""


VOTE_SYMBOLS = {
    Vote.UP: '<i class="fa fa-thumbs-o-up"></i>',
    Vote.BOOKMARK: '<i class="fa fa-bookmark"></i>',
    Vote.ACCEPT: '<i class="fa fa-heart"></i>',
}


@register.simple_tag
def vote_symbol(vote):
    return VOTE_SYMBOLS.get(vote.type, '')


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
    gravatar_url += urlencode(dict(s=str(size), d='identicon'))

    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, name)


# this contains the body of each comment
COMMENT_TEMPLATE = 'post_comment.html'
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
        # cont.update(csrf(request))
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


class MarkDownNode(template.Node):
    def __init__(self, nodelist):
        self.nodelist = nodelist


    def render(self, context):
        text = self.nodelist.render(context)
        text = html.sanitize(text, safe=True)
        return text


@register.tag('markdown')
def markdown_tag(parser, token):
    """
    Enables a block of markdown text to be used in a template.

    Syntax::

            {% markdown %}
            ## Markdown

            Now you can write markdown in your templates. This is good because:

            * markdown is awesome
            * markdown is less verbose than writing html by hand

            {% endmarkdown %}
    """
    nodelist = parser.parse(('endmarkdown',))
    # Need to do this otherwise we get big fail.
    parser.delete_first_token()
    return MarkDownNode(nodelist)