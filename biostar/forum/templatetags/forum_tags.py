import datetime
import itertools
import logging
import os
import random
from datetime import timedelta
from itertools import count, islice

import bleach
from django import template, forms
from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.paginator import Paginator
from django.db.models import Count
from django.shortcuts import reverse
from django.utils.safestring import mark_safe
from django.utils.timezone import utc
from taggit.models import Tag
from re import IGNORECASE, compile, escape
import html2markdown

from biostar.accounts.models import Profile, Message
from biostar.forum import const, auth
from biostar.utils import helpers
from biostar.forum import markdown
from biostar.forum.models import Post, Vote, Award, Subscription, Badge

User = get_user_model()

logger = logging.getLogger("engine")

register = template.Library()

ICON_MAP = dict(
    rank="list ol icon",
    views="eye icon",
    replies="comments icon",
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

@register.inclusion_tag('widgets/count_badge.html')
def count_badge(count):

    try:
        count = int(count)
    except ValueError as exc:
        # TODO: this is to catch ongoing stale sessions, may be removed later
        logger.info(f"invalid count (stale session?) {count}")
        count = 0

    return dict(count=count)


@register.simple_tag(takes_context=True)
def activate(context, state, target):
    targets = target.split(',')
    label = "active" if state in targets else ""

    return label


@register.simple_tag()
def vote_icon(vote):
    if vote.type == Vote.ACCEPT:
        return "check circle"
    elif vote.type == Vote.UP:
        return "thumbs up"
    elif vote.type == Vote.BOOKMARK:
        return "bookmark"

    return "help"


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

@register.filter
def htmltomarkdown(text):
    return helpers.htmltomarkdown(text)

@register.inclusion_tag('widgets/post_details.html', takes_context=True)
def post_details(context, post, user, avatar=True):
    return dict(post=post, user=user, avatar=avatar, context=context)


@register.simple_tag
def post_type_display(post_type):
    mapper = dict(Post.TYPE_CHOICES)
    return mapper.get(post_type)


def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)


@register.simple_tag
def gravatar(user=None, user_uid=None, size=80):
    if user_uid and hasattr(user, 'profile'):
        user = User.objects.filter(profile__uid=user_uid).first()

    return auth.gravatar(user=user, size=size)


@register.filter
def embed(text):
    return markdown.parse(text, clean=True, escape=True)


@register.inclusion_tag('widgets/filter_dropdown.html', takes_context=True)
def filter_dropdown(context):

    request = context['request']
    order = request.GET.get("order", 'rank') or 'rank'
    limit = request.GET.get("limit", 'all') or 'all'
    # TODO will be refactored out.
    context.update(dict(order=order, limit=limit))
    return context


@register.inclusion_tag('widgets/user_icon.html', takes_context=True)
def user_icon(context, user=None, is_moderator=False, is_spammer=False, score=0):
    try:
        is_moderator = user.profile.is_moderator if user else is_moderator
        score = user.profile.get_score() if user else score * 10
        is_spammer = user.profile.is_spammer if user else is_spammer
    except Exception as exc:
        logger.info(exc)

    context.update(dict(is_moderator=is_moderator, is_spammer=is_spammer, score=score))
    return context


@register.simple_tag()
def user_icon_css(user=None):
    css = ''
    if user and user.is_authenticated:

        if user.profile.is_moderator:
            css = "bolt icon"
        elif user.profile.score > 1000:
            css = "user icon"
        else:
            css = "user outline icon"

    return css


@register.inclusion_tag('widgets/post_user_line.html', takes_context=True)
def post_user_line(context, post, avatar=False, user_info=True):
    context.update(dict(post=post, avatar=avatar, user_info=user_info))
    return context


@register.inclusion_tag('widgets/post_user_line.html', takes_context=True)
def postuid_user_line(context, uid, avatar=True, user_info=True):
    post = Post.objects.filter(uid=uid).first()

    context.update(dict(post=post, avatar=avatar, user_info=user_info))
    return context


@register.inclusion_tag('widgets/user_card.html', takes_context=True)
def user_card(context, target):
    context.update(dict(target=target))
    return context


@register.inclusion_tag('widgets/post_user_box.html', takes_context=True)
def post_user_box(context, target_user, post):
    context.update(dict(target_user=target_user, post=post))
    return context


@register.inclusion_tag('widgets/post_actions.html', takes_context=True)
def post_actions(context, post, label="ADD COMMENT", author=None, lastedit_user=None, avatar=False):
    request = context["request"]

    return dict(post=post, user=request.user, author=author, lastedit_user=lastedit_user,
                label=label, request=request, avatar=avatar)


@register.inclusion_tag('widgets/post_tags.html')
def post_tags(post=None, post_uid=None, show_views=False, tags_str='', spaced=True):
    if post_uid:
        post = Post.objects.filter(uid=post_uid).first()

    tags = tags_str.split(",") if tags_str else ''
    tags = post.tag_val.split(",") if post else tags

    return dict(post=post, tags=tags, show_views=show_views, spaced=spaced)


@register.inclusion_tag('widgets/pages.html', takes_context=True)
def pages(context, objs, show_step=True):
    request = context["request"]
    url = request.path

    return dict(objs=objs, url=url, show_step=show_step, request=request)


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


@register.inclusion_tag('forms/field_tags.html', takes_context=True)
def tags_field(context, form_field, initial=''):
    """Render multi-select dropdown options for tags. """

    # Read from tags file
    tags_file = getattr(settings, "TAGS_OPTIONS_FILE", None)
    # Get currently selected tags from the post or request
    selected = initial.split(",") if initial else []
    selected = {(val, True) for val in selected}
    if tags_file:
        opts = file_tags_options(selected)
    else:
        opts = {}

    options = itertools.chain(selected, opts)

    context = dict(initial=initial, form_field=form_field, dropdown_options=options)

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


def get_tags_file():
    """
    Get a list of files to render from a file
    """
    # Get the tags op

    tags_file = getattr(settings, "TAGS_OPTIONS_FILE", None)

    return tags_file


def read_tags(exclude=[], limit=5000):
    """Read tags from a file. Each line is considered a tag. """
    # Get tags from a file
    tags_file = get_tags_file() or ''
    stream = open(tags_file, 'r') if os.path.exists(tags_file) else []
    stream = islice(zip(count(1), stream), limit)
    tags_opts = set()

    for idx, line in stream:
        line = line.strip()
        if line in exclude:
            continue
        tags_opts.add((line, False))

    return tags_opts


def file_tags_options(selected):
    """
    Present tags tags in a multi-select dropdown format.
    """
    # Read tags from file.
    try:
        opts = read_tags(exclude=selected)
    except Exception as exc:
        logger.error(f"Error reading tags from file.:{exc}")
        opts = []

    return opts


@register.inclusion_tag('forms/form_errors.html')
def form_errors(form, wmd_prefix='', override_content=False):
    """
    Turns form errors into a data structure
    """

    try:
        errorlist = [('', message, '') for message in form.non_field_errors()]
        for field in form:
            for error in field.errors:
                # wmd_prefix is required when dealing with 'content' field.
                field_id = wmd_prefix if (override_content and field.name == 'content') else field.id_for_label
                errorlist.append((f'{field.name}:', error, field_id))

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
def get_user_location(user):
    return user.profile.location or "location unknown"


@register.filter
def get_last_login(user):
    if user.profile.last_login:
        return f"{time_ago(user.profile.last_login)}"
    return f"{time_ago(user.profile.date_joined)}"


@register.inclusion_tag('widgets/feed_custom.html')
def custom_feed(objs, ftype='', title=''):
    users = ()
    if ftype == 'messages':
        users = set(m.sender for m in objs)
    if ftype in ['follow', 'bookmark', 'votes']:
        users = set(o.author for o in objs)

    context = dict(users=users, title=title)
    return context


@register.inclusion_tag(takes_context=True, filename='search/search_pages.html')
def search_pages(context, results):
    previous_page = results.pagenum - 1
    next_page = results.pagenum + 1 if not results.is_last_page() else results.pagenum
    request = context['request']
    query = request.GET.get('query', '')
    order = request.GET.get('order', 'relevance')
    context = dict(results=results, previous_page=previous_page, query=query,next_page=next_page, order=order)
    return context


@register.inclusion_tag(takes_context=True, filename='search/search_bar.html')
def search_bar(context, tags=False, users=False):
    search_url = reverse('tags_list') if tags else reverse('community_list') if users else reverse('post_search')
    request = context['request']
    value = request.GET.get('query', '')
    context = dict(search_url=search_url, value=value)

    return context


@register.simple_tag
def get_post_list(target, request, show=None):
    """
    Return post list belonging to a user
    """
    user = request.user
    page = request.GET.get("page", 1)
    posts = Post.objects.valid_posts(u=user, author=target)

    # Show a specific post listing.
    show_map = dict(questions=Post.QUESTION, tools=Post.TOOL, news=Post.NEWS,
                    blogs=Post.BLOG, tutorials=Post.TUTORIAL, answers=Post.ANSWER,
                    comments=Post.COMMENT, forum=Post.FORUM)

    type_filter = show_map.get(show)
    posts = posts.filter(type=type_filter) if type_filter is not None else posts

    posts = posts.select_related("root").select_related("author__profile", "lastedit_user__profile")
    posts = posts.order_by("-rank")

    # Cache the users posts add pagination to posts.
    paginator = Paginator(object_list=posts, per_page=settings.POSTS_PER_PAGE)
    posts = paginator.get_page(page)

    return posts


def awards_feed():
    awards = Award.objects.order_by('-pk').select_related("badge", "user", "user__profile")[:300]
    # Store already seen users
    seen = set()
    # Aggregate list of shown awards
    recent = []
    for award in awards:
        if award.user not in seen:
            recent.append(award)

        seen.update([award.user])
    recent = recent[:settings.AWARDS_FEED_COUNT]

    return recent


@register.inclusion_tag('widgets/feed_default.html')
def default_feed(user):
    recent_votes = Vote.objects.filter(post__status=Post.OPEN,
                                       post__root__status=Post.OPEN).prefetch_related("post")
    recent_votes = recent_votes.order_by("-pk")[:settings.VOTE_FEED_COUNT]

    # Get valid users that have a location set in profile.
    recent_locations = Profile.objects.valid_users().exclude(location="").prefetch_related("user")
    recent_locations = recent_locations.order_by('-last_login')[:settings.LOCATION_FEED_COUNT]

    # Adding the awards amount
    recent_awards = awards_feed()

    # Get valid posts
    recent_replies = Post.objects.valid_posts(is_toplevel=False).select_related("author__profile", "author")
    recent_replies = recent_replies.order_by("-pk")[:settings.REPLIES_FEED_COUNT]

    context = dict(recent_votes=recent_votes, recent_awards=recent_awards,
                   recent_locations=recent_locations, recent_replies=recent_replies,
                   user=user)

    return context


@register.simple_tag
def planet_gravatar(planet_author):
    email = planet_author.replace(' ', '')
    email = f"{email}@planet.org"
    email = email.encode('utf-8')
    return auth.gravatar_url(email=email, style='retro')


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


@register.inclusion_tag('widgets/list_badges.html', takes_context=True)
def list_awards(context, target):
    request = context['request']
    # Show list of all awards here.
    badges = Badge.objects.filter(award__user=target).annotate(count=Count("award"))
    badges = badges.order_by('-count')
    context = dict(badges=badges, request=request, target=target)
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


@register.simple_tag(takes_context=True)
def relative_url(context, value, field_name, urlencode=None):
    """
    Updates field_name parameters in url with new value
    """
    # Check if the param is in biggest common set.
    expect = const.ALLOWED_PARAMS

    def apply_filter(param):
        key = param.split('=')[0]
        # Return parameter if in valid const.
        if key != field_name and key in expect:
            return param

    request = context['request']
    # Create more_like_this string with updated field_name, value pair.
    url = f'?{field_name}={value}'
    if urlencode:
        # Split query string
        querystring = urlencode.split('&')
        # Exclude old value 'field_name' from query string
        filtered = filter(apply_filter, querystring)
        # Join the filtered string
        encoded = '&'.join(filtered)
        # Update query string
        encoded = f'&{encoded}' if encoded else ''
        url = f'{url}{encoded}'

    # Concatenate current path to relative query string
    url = request.path + url

    return url


@register.simple_tag
def get_thread_users(users, post, limit=2):
    displayed_users = {post.author, post.lastedit_user or post.author}

    for user in users:
        if len(displayed_users) >= limit:
            break
        if user in displayed_users:
            continue

        displayed_users.add(user)

    return displayed_users


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


def post_boxclass(root_type, answer_count, root_has_accepted):
    # Create the css class for each row
    if root_type == Post.JOB:
        style = "job"
    elif root_type == Post.TUTORIAL:
        style = "tutorial"
    elif root_type == Post.TOOL:
        style = "tool"
    elif root_type == Post.FORUM:
        style = "forum"
    elif root_type == Post.NEWS:
        style = "news"
    elif root_type == Post.HERALD:
        style = "herald"
    else:
        style = "question"

    if isinstance(answer_count, int) and int(answer_count) >= 1:
        style += " has_answers"

    if root_has_accepted == True:
        modifier = "accepted answer" if root_type == Post.QUESTION else "accepted"
    else:
        modifier = "open"

    return f"{style} {modifier}"


@register.inclusion_tag('herald/herald_item.html', takes_context=True)
def herald_item(context, item):
    request = context['request']
    user = request.user
    context = dict(story=item, request=request, user=user)
    return context


@register.simple_tag
def search_boxclass(root_type, answer_count, root_has_accepted):
    return post_boxclass(root_type=root_type,
                         answer_count=answer_count,
                         root_has_accepted=root_has_accepted)


@register.simple_tag
def boxclass(post=None, uid=None):
    return post_boxclass(root_type=post.root.type,
                         answer_count=post.root.answer_count,
                         root_has_accepted=post.root.has_accepted)


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
        collect.append(f'<div class="indent" ><div>{html}</div>')

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


def top_level_only(attrs, new=False):
    '''
    Helper function used when linkifying with bleach.
    '''
    if not new:
        return attrs
    text = attrs['_text']
    if not text.startswith(('http:', 'https:')):
        return None
    return attrs


@register.simple_tag
def markdown_file(pattern):
    """
    Returns the content of a file matched by the pattern.
    Returns an error message if the pattern cannot be found.
    """
    path = pattern
    path = os.path.abspath(path)
    if os.path.isfile(path):
        text = open(path).read()
    else:
        text = f" file '{pattern}': '{path}' not found"

    try:

        html = markdown.parse(text, clean=False, escape=False, allow_rewrite=True)
        html = bleach.linkify(html, callbacks=[top_level_only], skip_tags=['pre'])
        html = mark_safe(html)
    except Exception as e:
        html = f"Markdown rendering exception"
        logger.error(e)
    return html


class MarkDownNode(template.Node):
    # CALLBACKS = [top_level_only]

    def __init__(self, nodelist):
        self.nodelist = nodelist

    def render(self, context):
        text = self.nodelist.render(context)
        text = markdown.parse(text, clean=False, escape=False, allow_rewrite=True)
        # text = bleach.linkify(text, callbacks=self.CALLBACKS, skip_tags=['pre'])
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

    # need to do this otherwise we get big fail
    parser.delete_first_token()
    return MarkDownNode(nodelist)
