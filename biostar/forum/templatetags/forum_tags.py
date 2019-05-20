import hashlib
import itertools
import logging
import random
import urllib.parse
from datetime import datetime
from datetime import timedelta

from django import template
from django.db.models import Q
from django.conf import settings
from django.contrib.auth import get_user_model
from django.utils.safestring import mark_safe
from django.utils.timezone import utc

from biostar.forum import forms, models, const, util
from biostar.forum.models import Post, Vote, Award
from biostar.message.models import Message
from biostar.utils.shortcuts import reverse

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


@register.simple_tag
def activate(state, target):
    return "active" if state == target else ""


def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.inclusion_tag('widgets/user_box.html')
def user_box(user, post):
    return dict(user=user, post=post)


@register.simple_tag
def get_all_message_count(request):
    user = request.user
    outbox = inbox = projects = mentioned = unread = 0

    if user.is_authenticated:
        inbox = Message.objects.inbox_for(user=user)
        outbox = Message.objects.outbox_for(user=user).count()
        unread = inbox.filter(unread=True).count()
        mentioned = inbox.filter(source=Message.MENTIONED).count()

    context = dict(outbox=outbox, inbox=inbox.count(), projects=projects, mentioned=mentioned,
                   unread=unread)

    return context


@register.inclusion_tag('widgets/pages.html')
def pages(objs, request):
    url = request.path

    return dict(objs=objs, url=url, request=request)


@register.simple_tag
def get_tags_list(tags_str):
    return set(util.split_tags(tags_str))


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


@register.simple_tag
def gravatar(user, size=80):
    # name = user.profile.name
    if user.is_anonymous or user.profile.is_suspended:
        # Removes spammy images for suspended users
        email = 'suspended@biostars.org'.encode('utf8')
    else:
        email = user.email.encode('utf8')

    hash = hashlib.md5(email).hexdigest()

    url = "https://secure.gravatar.com/avatar/%s?" % hash
    url += urllib.parse.urlencode({
        's': str(size),
        'd': 'retro',
    })
    return url


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree, form):
    "Renders the post body"
    request = context['request']
    return dict(post=post, user=user, tree=tree, request=request,
                form=form,
                redir_field_name=const.REDIRECT_FIELD_NAME)


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


@register.filter
def show_email(user):
    try:
        head, tail = user.email.split("@")
        email = head[0] + "*" * 10 + tail
    except:
        return user.email[0] + "*" * 10

    return email


@register.simple_tag
def is_moderator(user):
    if user.is_authenticated and user.profile.is_moderator:
        return True
    return False


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

    posts = Post.objects.exclude(uid=post.uid).filter(query)[:settings.SINGLE_FEED_COUNT]
    context = dict(posts=posts)
    return context


@register.inclusion_tag('widgets/feed.html')
def feed(user):
    recent_votes = Vote.objects.prefetch_related("post")[:settings.VOTE_FEED_COUNT]

    recent_locations = User.objects.exclude(profile__location="")
    recent_locations = recent_locations.prefetch_related("profile")[:settings.LOCATION_FEED_COUNT]

    recent_awards = Award.objects.order_by("-pk").select_related("badge", "user", "user__profile")
    recent_awards = recent_awards[:settings.AWARDS_FEED_COUNT]
    recent_replies = Post.objects.filter(type__in=[Post.COMMENT, Post.ANSWER])
    recent_replies = recent_replies.select_related("author__profile", "author")[:settings.REPLIES_FEED_COUNT]

    context = dict(recent_votes=recent_votes, recent_awards=recent_awards,
                   recent_locations=recent_locations, recent_replies=recent_replies,
                   user=user)

    return context


@register.filter
def show_score_icon(user):
    return icon(user)


@register.filter
def icon(user):
    if user.profile.is_moderator:
        icon = f'<i class="ui muted bolt icon"></i>'
    elif user.profile.score > 100:
        icon = f'<i class="ui muted user icon"></i>'
    else:
        icon = f'<i class="ui muted user outline icon"></i>'

    return mark_safe(icon)


@register.filter
def show_score(score):
    score = (score * 10) + 1
    return score


@register.inclusion_tag('widgets/render_tags.html')
def render_tags(post):
    tags = post.tag_val.split(",")
    return dict(tags=tags)


@register.inclusion_tag('widgets/user_info.html')
def user_info(post, by_diff=False, with_image=True):
    return dict(post=post, by_diff=by_diff, with_image=with_image)


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


@register.inclusion_tag('widgets/listing.html')
def listing(post=None, user=None):
    if user is None:
        return dict(post=post)

    if post.is_deleted and (user.is_anonymous or not user.profile.is_moderator):
        return None

    return dict(post=post)


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
    elif post.has_accepted:
        style = "accept"
    elif post.reply_count > 0:
        style = "lightgreen"
    elif post.comment_count > 0:
        style = "grey"
    else:
        style = "maroon"

    return style


@register.simple_tag(takes_context=True)
def render_comments(context, tree, post, comment_template='widgets/comment_body.html'):
    request = context["request"]
    if post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree, comment_template=comment_template)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, comment_template):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(comment_template)

    def traverse(node):

        data = ['<div class="ui comment segments">']
        cont = {"post": node, 'user': request.user, 'request': request}
        html = body.render(cont)
        data.append(html)
        for child in tree.get(node.id, []):
            data.append(f'<div class="ui segment comment basic">')
            data.append(traverse(child))
            data.append("</div>")

        data.append("</div>")
        return '\n'.join(data)

    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:
        coll.append(traverse(node))

    return '\n'.join(coll)
