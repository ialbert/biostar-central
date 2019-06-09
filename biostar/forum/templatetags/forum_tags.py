import hashlib
import itertools
import logging
import random
from datetime import datetime
from datetime import timedelta

import hashlib
import urllib.parse

from django import template
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db.models import Q
from django.utils.safestring import mark_safe
from django.utils.timezone import utc

from biostar.accounts.models import Profile
from biostar.forum import const
from biostar.forum.models import Post, Vote, Award, Subscription

User = get_user_model()

logger = logging.getLogger("engine")

register = template.Library()

ICON_MAP = dict(
    rank="list ol icon",
    views="eye icon",
    replies="comment icon",
    votes="thumbs up icon",
    all='calendar plus icon',
    today='clock icon',
    week='calendar minus outline icon',
    month='calendar alternate icon',
    year='calendar icon',
    visit='sort numeric down icon',
    reputation='star icon',
    joined='sign in icon',
    activity='comment icon',
    rsent="sort numeric down icon",
    sent="sort numeric up icon",
    rep="user outline icon"
)


@register.simple_tag(takes_context=True)
def activate(context, state, target):
    label = "active" if state == target else ""
    request = context['request']
    count = 0

    # Special casing a few targets to generate an extra css class.
    if target == "messages":
        count = request.session.get("counts", {}).get("message_count", 0)
    elif target == "votes":
        count = request.session.get("counts", {}).get("vote_count", 0)

    # Generate a broader css if necessary.
    label = f"new {label}" if count else label

    return label

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


@register.simple_tag(takes_context=True)
def count_label(context, label):
    request = context['request']
    count = request.session.get("counts", {}).get(label, 0)
    label = f"({count})" if count else ""
    return label


@register.inclusion_tag('widgets/inplace_form.html')
def inplace_form(post, width='100%'):

    pad = 4 if post.type == Post.COMMENT else 7
    rows = len(post.content.split("\n")) + pad
    context = dict(post=post, width=width, rows=rows)
    return context


def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.simple_tag
def gravatar(user, size=80):
    style = "retro"
    if user.is_anonymous or user.profile.is_suspended:
        # Removes spammy images for suspended users
        # email = 'suspended@biostars.org'.encode('utf8')
        style = "monsterid"
    else:
        if user.profile.is_moderator:
            style = "robohash"
        email = user.email.encode('utf8')

    hash = hashlib.md5(email).hexdigest()

    gravatar_url = "https://secure.gravatar.com/avatar/%s?" % hash
    gravatar_url += urllib.parse.urlencode({
        's': str(size),
        'd': style,
    }
    )
    return gravatar_url


@register.simple_tag()
def user_score(user):
    score = user.profile.score * 10
    return score


@register.inclusion_tag('widgets/user_icon.html')
def user_icon(user):
    score = user_score(user)
    context = dict(user=user, score=score)
    return context

@register.inclusion_tag('widgets/post_user_line.html')
def post_user_line(post, avatar=False):
    return dict(post=post, avatar=avatar)


@register.inclusion_tag('widgets/user_card.html')
def user_card(user):
    return dict(user=user)


@register.inclusion_tag('widgets/post_user_box.html')
def post_user_box(user, post):
    return dict(user=user, post=post)


@register.inclusion_tag('widgets/post_actions.html', takes_context=True)
def post_actions(context, post, label="ADD COMMENT", avatar=False):
    request = context["request"]
    return dict(post=post, label=label, request=request, avatar=avatar)


@register.inclusion_tag('widgets/post_tags.html')
def post_tags(post, show_views=False):
    tags = post.tag_val.split(",")
    return dict(post=post, tags=tags, show_views=show_views)


@register.inclusion_tag('widgets/pages.html', takes_context=True)
def pages(context, objs):
    request = context["request"]
    url = request.path
    return dict(objs=objs, url=url, request=request)


@register.simple_tag
def randparam():
    "Append to URL to bypass server caching of CSS or JS files"
    return f"?randval={random.randint(1, 10000000)}" if settings.DEBUG else ""


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.filter
def unread(message, user):
    if message.recipient == user and message.unread:
        return "unread-message"
    return ""


@register.simple_tag(takes_context=True)
def follow_label(context, post):
    user = context["request"].user

    not_following = "not following"

    label_map = {
        Profile.LOCAL_MESSAGE: "following with messages",
        Profile.DEFAULT_MESSAGES: "following with messages",
        Profile.EMAIL_MESSAGE: "following via email",
        Profile.NO_MESSAGES: not_following
    }

    if user.is_anonymous:
        return not_following

    # Get the current subscription
    sub = Subscription.objects.filter(post=post, user=user).first()
    sub = sub or Subscription(post=post, user=user, type=Profile.NO_MESSAGES)

    label = label_map.get(sub.type, not_following)

    return label


@register.inclusion_tag('widgets/form_errors.html')
def form_errors(form):
    """
    Turns form errors into a data structure
    """

    try:
        errorlist = [('', message) for message in form.non_field_errors()]
        for field in form:
            for error in field.errors:
                errorlist.append((f'{field.name}:', error, field.id_for_label))
    except Exception as exc:
        errorlist = []
        logging.error(exc)

    context = dict(errorlist=errorlist)

    return context


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree):
    "Renders the post body"
    request = context['request']
    return dict(post=post, user=user, tree=tree, request=request)


@register.filter
def get_award_context(award):
    post = award.post
    context = f"For : <a href={post.get_absolute_url()}>{post.title}</a>" if post else ""
    return context


@register.filter
def get_user_location(user):
    return user.profile.location or "location unknown"


@register.filter
def get_last_login(user):
    if user.profile.last_login:
        return f"{time_ago(user.profile.last_login)}"
    return f"{time_ago(user.profile.date_joined)}"


@register.inclusion_tag('widgets/single_feed.html')
def single_post_feed(post):
    """
    Return single post feed populated with similar posts.
    """
    tags = post.tag_val.split(",")

    # Gather similar posts
    query = Q()
    for tag in tags:
        query |= Q(tag_val__iregex=tag)

    posts = Post.objects.exclude(uid=post.uid).filter(query, type__in=Post.TOP_LEVEL)[:settings.SINGLE_FEED_COUNT]
    context = dict(posts=posts)
    return context


@register.inclusion_tag('widgets/listing.html', takes_context=True)
def list_posts(context, user):
    posts = Post.objects.filter(author=user)
    request = context["request"]
    context = dict(posts=posts, request=request)
    return context


@register.inclusion_tag('widgets/feed.html')
def feed(user):
    recent_votes = Vote.objects.prefetch_related("post")
    recent_votes = recent_votes.order_by("-pk")[:settings.VOTE_FEED_COUNT]

    recent_locations = User.objects.exclude(profile__location="")
    recent_locations = recent_locations.prefetch_related("profile")[:settings.LOCATION_FEED_COUNT]

    recent_awards = Award.objects.order_by("-pk").select_related("badge", "user", "user__profile")
    recent_awards = recent_awards[:settings.AWARDS_FEED_COUNT]

    recent_replies = Post.objects.filter(type__in=[Post.COMMENT, Post.ANSWER])
    recent_replies = recent_replies.select_related("author__profile", "author")
    recent_replies = recent_replies.order_by("-pk")[:settings.REPLIES_FEED_COUNT]

    context = dict(recent_votes=recent_votes, recent_awards=recent_awards,
                   recent_locations=recent_locations, recent_replies=recent_replies,
                   user=user)

    return context


@register.simple_tag
def get_icon(string, default=""):
    icon = ICON_MAP.get(string) or ICON_MAP.get(default)
    return icon


@register.simple_tag
def get_wording(filtered, prefix="Sort by:", default=""):
    """
    Get the naming and icons for limits and ordering.
    """

    display = dict(all="all time", week="this week", month="this month",
                   year="this year", rank="rank", views="views", today="today",
                   replies="replies", votes="votes", visit="recent visit",
                   reputation="reputation", joined="date joined", activity="activity level",
                   rsent="oldest to newest ", sent="newest to oldest",
                   rep="sender reputation")
    if display.get(filtered):
        displayed = display[filtered]
    else:
        displayed = display[default]

    wording = f"{prefix} {displayed}"

    return wording


@register.simple_tag
def relative_url(value, field_name, urlencode=None):
    """
    Updates field_name parameters in url with value
    """
    # Create query string with updated field_name, value pair.
    url = '?{}={}'.format(field_name, value)
    if urlencode:
        # Split query string
        querystring = urlencode.split('&')
        # Exclude old value 'field_name' from query string
        filter_func = lambda p: p.split('=')[0] != field_name
        filtered_querystring = filter(filter_func, querystring)
        # Join the filtered string
        encoded_querystring = '&'.join(filtered_querystring)
        # Update query string
        url = '{}&{}'.format(url, encoded_querystring)

    return url


@register.simple_tag
def get_thread_users(post, limit=5):
    thread_users = post.thread_users.all()
    stream = itertools.islice(thread_users, limit)

    # Author is shown first
    users = [post.author]
    for user in stream:
        if user in users:
            continue
        users.append(user)

    return users


@register.inclusion_tag('widgets/listing.html', takes_context=True)
def listing(context, posts=None):
    request = context["request"]
    return dict(posts=posts, request=request)


@register.filter
def show_nonzero(value):
    "The purpose of this is to return value or empty"
    return value if value else ''


@register.simple_tag
def object_count(request, otype):
    user = request.user
    count = 0

    if user.is_authenticated:

        if otype == "message":
            count = user.profile.new_messages

    return count


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


@register.filter
def time_ago(date):
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
def get_subcount(sub_count):
    if sub_count > 5:
        return mark_safe(f'<div class="subs">{sub_count} follow</div>')

    return ""


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

    if post.type == Post.JOB:
        style = "job"
    elif post.type == Post.TUTORIAL:
        style = "tutorial"
    elif post.type == Post.TOOL:
        style = "tool"
    elif post.type == Post.FORUM:
        style = "forum"
    elif post.type == Post.NEWS:
        style = "news"
    elif post.answer_count:
        style = "has_answers"
    else:
        style = "question"

    if post.has_accepted:
        modifier = "accepted answer" if post.type == Post.QUESTION else "accepted"
    else:
        modifier = "open"

    return f"{style} {modifier}"


@register.simple_tag(takes_context=True)
def render_comments(context, tree, post, template_name='widgets/comment_body.html'):
    request = context["request"]
    if post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree, template_name=template_name)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, template_name):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(template_name)

    seen = set()

    def traverse(node, collect=[]):

        cont = {"post": node, 'user': request.user, 'request': request}
        html = body.render(cont)

        collect.append(f'<div class="indent"><div class="comment">{html}</div>')

        for child in tree.get(node.id, []):
            if child in seen:
                raise Exception(f"circular tree {child.pk} {child.title}")
            seen.add(child)
            traverse(child, collect=collect)

        collect.append(f"</div>")

    # this collects the comments for the post
    collect = ['<div class="comment-list">']
    for node in tree[post.id]:
        traverse(node, collect=collect)
    collect.append("</div>")

    html = '\n'.join(collect)

    return html
