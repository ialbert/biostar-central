import hashlib
import itertools
import logging
import random
import re
import urllib.parse
import datetime
from itertools import count, islice
from datetime import timedelta

from snowpenguin.django.recaptcha2.fields import ReCaptchaField
from snowpenguin.django.recaptcha2.widgets import ReCaptchaWidget
from taggit.models import Tag
from django import template, forms
from django.db.models import Count
from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.paginator import Paginator
from django.shortcuts import reverse
from django.db.models import Q
from django.utils.safestring import mark_safe
from django.utils.timezone import utc

from biostar.accounts.models import Profile, Message
from biostar.forum import const, auth
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
    rep="user outline icon",
    tagged="tags icon",
)


def get_count(request, key, default=0):
    """
    Returns a count stored in the session.
    """
    value = request.session.get(const.COUNT_DATA_KEY, {}).get(key, default)
    return value


@register.simple_tag(takes_context=True)
def activate(context, state, target):
    label = "active" if state == target else ""
    request = context['request']
    value = 0

    # Special casing a few targets to generate an extra css class.
    if target == "messages":
        value = get_count(request, "message_count")
    elif target == "votes":
        value = get_count(request, "vote_count")

    # Generate a broader css if necessary.
    label = f"new {label}" if value else label

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
    count = get_count(request, label)
    label = f"({count})" if count else ""
    return label


@register.simple_tag
def users_list():
    users = User.objects.exclude(profile__state__in=[Profile.DEACTIVATED, Profile.SUSPENDED, Profile.BANNED])
    users = ','.join([x.username for x in users])

    return users


@register.inclusion_tag('widgets/inplace_form.html')
def inplace_form(post, width='100%'):
    pad = 4 if post.type == Post.COMMENT else 7
    rows = len(post.content.split("\n")) + pad
    context = dict(post=post, width=width, rows=rows)
    return context


@register.inclusion_tag('widgets/post_user_line.html')
def post_search_line(post_uid, avatar=True):
    post = Post.objects.filter(uid=post_uid).first()
    return dict(post=post, avatar=avatar)


@register.inclusion_tag('widgets/pages_search.html', takes_context=True)
def pages_search(context, results):

    previous_page = results.pagenum - 1
    next_page = results.pagenum + 1 if not results.is_last_page() else results.pagenum
    request = context['request']
    query = request.GET.get('query', '')
    context = dict(results=results, previous_page=previous_page, query=query,
                   next_page=next_page)

    return context


def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)


@register.simple_tag
def gravatar(user, size=80):

    return  auth.gravatar(user=user, size=size)

@register.simple_tag()
def user_score(user):
    score = user.profile.score * 10
    return score


@register.inclusion_tag('widgets/user_icon.html')
def user_icon(user=None, user_uid=None):
    try:
        user = user or User.objects.filter(profile__uid=user_uid).first()
        score = user_score(user)
    except Exception as exc:
        logger.info(exc)
        user = score = None

    context = dict(user=user, score=score)
    return context


@register.inclusion_tag('widgets/post_user_line.html')
def post_user_line(post, avatar=False, user_info=True):
    return dict(post=post, avatar=avatar, user_info=user_info)

@register.inclusion_tag('widgets/post_user_line.html')
def postuid_user_line(uid, avatar=True, user_info=True):
    post = Post.objects.filter(uid=uid).first()
    return dict(post=post, avatar=avatar, user_info=user_info)


@register.inclusion_tag('widgets/user_card.html')
def user_card(user):
    return dict(user=user)


@register.inclusion_tag('widgets/post_user_box.html')
def post_user_box(user, post):
    return dict(user=user, post=post)


@register.inclusion_tag('widgets/post_actions.html', takes_context=True)
def post_actions(context, post, label="ADD COMMENT", avatar=False):
    request = context["request"]

    return dict(post=post, user=request.user,
                label=label, request=request, avatar=avatar)


@register.inclusion_tag('widgets/post_tags.html')
def post_tags(post=None, post_uid=None, show_views=False, spaced=True):

    if post_uid:
        post = Post.objects.filter(uid=post_uid).first()

    tags = post.tag_val.split(",")

    return dict(post=post, tags=tags, show_views=show_views, spaced=spaced)


@register.simple_tag
def get_vote_count(uid):
    post = Post.objects.filter(uid=uid).first()
    return post.get_votecount


@register.simple_tag
def get_view_count(uid):
    post = Post.objects.filter(uid=uid).first()
    return post.root.view_count


@register.simple_tag
def get_subs_count(uid):
    post = Post.objects.filter(uid=uid).first()
    return post.subs_count


@register.simple_tag
def get_reply_count(uid):
    post = Post.objects.filter(uid=uid).first()
    return post.reply_count


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


@register.simple_tag
def toggle_unread(user):
    Message.objects.filter(recipient=user, unread=True).update(unread=False)
    return ''


@register.simple_tag(takes_context=True)
def digest_label(context, post):

    user = context['request'].user
    no_digest = 'No digest'

    label_map = {
        Profile.WEEKLY_DIGEST: "Weekly digest",
        Profile.MONTHLY_DIGEST: "Monthly digest",
        Profile.DAILY_DIGEST: 'Daily digest',
        Profile.NO_DIGEST: no_digest
    }
    if user.is_anonymous:
        return no_digest

    label = label_map.get(user.profile.digest_prefs, no_digest)

    return label


@register.simple_tag(takes_context=True)
def follow_label(context, post):
    user = context["request"].user

    not_following = "not following"

    label_map = {
        Subscription.LOCAL_MESSAGE: "following with messages",
        Subscription.EMAIL_MESSAGE: "following via email",
        Subscription.NO_MESSAGES: not_following,
    }

    if user.is_anonymous:
        return not_following

    # Get the current subscription
    sub = Subscription.objects.filter(post=post.root, user=user).first()
    sub = sub or Subscription(post=post, user=user, type=Subscription.NO_MESSAGES)

    label = label_map.get(sub.type, not_following)

    return label


@register.inclusion_tag('widgets/type_help.html')
def type_help():

    context = dict()
    return context

@register.simple_tag
def inplace_type_field(post=None, field_id='type'):
    choices = [opt for opt in Post.TYPE_CHOICES]

    choices = filter(lambda opt: (opt[1] in settings.ALLOWED_POST_TYPES) if settings.ALLOWED_POST_TYPES else
                                 (opt[0] in Post.TOP_LEVEL), choices)

    post_type = forms.IntegerField(label="Post Type",
                                   widget=forms.Select(choices=choices, attrs={'class': "ui fluid dropdown",
                                                                               'id': field_id}),
                                   help_text="Select a post type.")

    value = post.type if post else Post.QUESTION
    post_type = post_type.widget.render('post_type', value)

    return mark_safe(post_type)


@register.simple_tag
def get_tags(request=None, post=None, user=None, watched=False):
    # Get tags in requests before fetching ones in the post.
    # This is done to accommodate populating tags in forms

    if request:
        tags = request.GET.get('tag_val', request.POST.get('tag_val', ''))
    else:
        tags = post.tag_val if isinstance(post, Post) else ''
        tags = user.profile.my_tags if user else tags
        tags = user.profile.watched_tags if watched else tags

    # Prepare the tags options in the dropdown from a file
    if settings.TAGS_OPTIONS_FILE:
        tags_opts = open(settings.TAGS_OPTIONS_FILE, 'r').readlines()
        tags_opts = [(x.strip(), False) if x.strip() not in tags.split(",") else (x.strip(), True)
                     for x in tags_opts if x != '\n']
    # Prepare dropdown options from database.
    else:
        query = Count('post')
        tags_query = Tag.objects.annotate(count=query).order_by('-count').exclude(name__in=tags.split(','))[:50]
        tags_opts = ((tag.name.strip(), False) for tag in tags_query)

    selected_tags_opt = ((val, True) for val in tags.split(","))
    tags_opts = itertools.chain(selected_tags_opt, tags_opts)

    context = dict(selected=tags, tags_opt=tags_opts)

    return context


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


@register.filter
def highlight(hit, field):
    lit = hit.highlights(field, top=5)
    return mark_safe(lit) if len(lit) else hit[field]


@register.inclusion_tag('widgets/feed_custom.html')
def custom_feed(objs, feed_type='', title=''):
    users = ()
    if feed_type == 'messages':
        users = set(m.sender for m in objs)
    if feed_type in ['following', 'bookmarks', 'votes']:
        users = set(o.author for o in objs)

    context = dict(users=users, title=title)
    return context


@register.inclusion_tag(takes_context=True, filename='widgets/search_bar.html')
def search_bar(context, search_url='', tags=False, users=False, ajax_results=True, extra_css='',
               redir=False):
    search_url = search_url or reverse('ajax_search')
    redir = '1' if redir else '0'
    request = context['request']
    value = request.GET.get('query', '')
    context = dict(search_url=search_url, tags=tags, users=users, extra_css=extra_css,
                   ajax_results=ajax_results, redir=redir, value=value)

    return context


@register.inclusion_tag('widgets/listing.html', takes_context=True)
def list_posts(context, target):
    request = context["request"]
    user = request.user
    posts = Post.objects.filter(author=target)
    page = request.GET.get('page', 1)
    posts = posts.prefetch_related("root", "author__profile",
                                   "lastedit_user__profile", "thread_users__profile")
    # Filter deleted items for anonymous and non-moderators.
    if user.is_anonymous or (user.is_authenticated and not user.profile.is_moderator):
        posts = posts.exclude(Q(status=Post.DELETED))
    posts = posts.order_by("-rank")
    posts = posts.exclude(Q(root=None) | Q(parent=None))
    # Create the paginator and apply post paging
    paginator = Paginator(posts, settings.POSTS_PER_PAGE)
    posts = paginator.get_page(page)

    request = context["request"]
    context = dict(posts=posts, request=request, include_pages_bar=True)
    return context


@register.inclusion_tag('widgets/feed_default.html')
def default_feed(user):
    recent_votes = Vote.objects.prefetch_related("post").exclude(post__status=Post.DELETED)
    recent_votes = recent_votes.order_by("-pk")[:settings.VOTE_FEED_COUNT]

    recent_locations = Profile.objects.exclude(Q(location="") | Q(state__in=[Profile.BANNED, Profile.SUSPENDED]))
    recent_locations = recent_locations.order_by('-last_login')
    recent_locations = recent_locations[:settings.LOCATION_FEED_COUNT]

    recent_awards = Award.objects.order_by("-pk").select_related("badge", "user", "user__profile")
    recent_awards = recent_awards.exclude(user__profile__state__in=[Profile.BANNED, Profile.SUSPENDED])
    recent_awards = recent_awards[:settings.AWARDS_FEED_COUNT]

    recent_replies = Post.objects.filter(type__in=[Post.COMMENT, Post.ANSWER]).exclude(status=Post.DELETED)
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
def get_digest_icon(user):
    no_digest = 'bell slash icon'

    icon_map = {Profile.WEEKLY_DIGEST: 'hourglass icon', Profile.MONTHLY_DIGEST: 'calendar icon',
                Profile.DAILY_DIGEST: 'clock icon', Profile.NO_DIGEST: no_digest}

    icon = icon_map.get(user.profile.digest_prefs) or no_digest
    return icon


@register.inclusion_tag('widgets/list_awards.html', takes_context=True)
def list_awards(context, target):
    request = context['request']
    awards = Award.objects.filter(user=target).order_by("-date")
    page = request.GET.get('page', 1)
    # Create the paginator
    paginator = Paginator(awards, 20)

    # Apply the votes paging.
    awards = paginator.get_page(page)

    context = dict(awards=awards, request=request)
    return context


@register.simple_tag
def get_wording(filtered, prefix="Sort by:", default=""):
    """
    Get the naming and icons for limits and ordering.
    """

    display = dict(all="all time", week="this week", month="this month",
                   year="this year", rank="Rank", views="Views", today="today",
                   replies="replies", votes="Votes", visit="recent visit",
                   reputation="reputation", joined="date joined", activity="activity level",
                   rsent="oldest to newest ", sent="newest to oldest",
                   rep="sender reputation", tagged="tagged")
    if display.get(filtered):
        displayed = display[filtered]
    else:
        displayed = display[default]

    wording = f"{prefix} {displayed}"

    return wording


@register.simple_tag
def activate_check_mark(filter, active):

    if filter == active:
        return 'check icon'

    return ''


@register.simple_tag
def relative_url(value, field_name, urlencode=None):
    """
    Updates field_name parameters in url with new value
    """
    # Create preform_search string with updated field_name, value pair.
    url = f'?{field_name}={value}'
    if urlencode:
        # Split preform_search string
        querystring = urlencode.split('&')
        # Exclude old value 'field_name' from preform_search string
        filter_func = lambda p: p.split('=')[0] != field_name
        filtered_querystring = filter(filter_func, querystring)
        # Join the filtered string
        encoded_querystring = '&'.join(filtered_querystring)
        # Update preform_search string
        url = f'{url}&{encoded_querystring}'

    return url


@register.simple_tag
def get_thread_users(post, limit=2):
    users = post.thread_users.exclude(profile__state__in=[Profile.BANNED, Profile.SUSPENDED])

    displayed_users = {post.author, post.lastedit_user}

    displayed_users = [u for u in displayed_users
                       if u.profile.state not in (Profile.BANNED, Profile.SUSPENDED)]
    for user in users:
        if len(displayed_users) >= limit:
            break
        if user in displayed_users:
            continue
        displayed_users.append(user)

    return displayed_users


@register.inclusion_tag('widgets/listing.html', takes_context=True)
def listing(context, posts=None, show_subs=True):
    request = context["request"]
    return dict(posts=posts, request=request, show_subs=show_subs)


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


@register.simple_tag
def subscription_label(sub_count):
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
def boxclass(post=None, uid=None):

    if uid:
        post = Post.objects.filter(uid=uid).first()

    # Create the css class for each row
    if post.root.type == Post.JOB:
        style = "job"
    elif post.root.type == Post.TUTORIAL:
        style = "tutorial"
    elif post.root.type == Post.TOOL:
        style = "tool"
    elif post.root.type == Post.FORUM:
        style = "forum"
    elif post.root.type == Post.NEWS:
        style = "news"
    else:
        style = "question"

    if post.root.answer_count:
        style += " has_answers"

    if post.root.has_accepted:
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
        source = f"indent-{node.uid}"
        target = f"'{node.uid}'"
        if request.user.is_authenticated and request.user.profile.is_moderator:
            collect.append(f'<div class="indent " id="{source}" ondragover="allowDrop(event);" ondrop="drop(event, {target})"><div class="comment">{html}</div>')
        else:
            collect.append(f'<div class="indent "><div class="comment">{html}</div>')

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


@register.filter
def get_children_list(post):

    children = []
    auth.walk_down_thread(parent=post, collect=children, is_root=post.is_toplevel)

    # Include itself in list

    children = [post.uid] + list(map(lambda p: p.uid, children))

    children += ['NEW'] if post.is_answer else []
    #print(children, post.uid)
    return ','.join(children)
