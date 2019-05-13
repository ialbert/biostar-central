import logging
import hashlib
import urllib.parse
import itertools

from datetime import timedelta, datetime
from django.utils.timezone import utc
from django import template
from django.utils.safestring import mark_safe
from django.core.paginator import Paginator
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db.models import Count
from django.db.models import Q
from datetime import datetime
from django.contrib import messages

from biostar.utils.shortcuts import reverse
from biostar.forum.models import Post, Vote, Award
from biostar.forum import auth, forms, models, const, util


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

def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.inclusion_tag('widgets/user_box.html')
def user_box(user):

    return dict(user=user)


@register.inclusion_tag('widgets/pages.html')
def pages(objs, request):

    url = request.path

    return dict(objs=objs, url=url, request=request)


@register.simple_tag
def get_tags_list(tags_str):

    return set(util.split_tags(tags_str))


@register.simple_tag
def gravatar(user, size=80):
    #name = user.profile.name
    if user.is_anonymous or user.profile.is_suspended:
        # Removes spammy images for suspended users
        email = 'suspended@biostars.org'.encode('utf8')
    else:
        email = user.email.encode('utf8')

    hash = hashlib.md5(email).hexdigest()

    gravatar_url = "https://secure.gravatar.com/avatar/%s?" % hash
    gravatar_url += urllib.parse.urlencode({
        's': str(size),
        'd': 'retro',
    }
    )

    return mark_safe(f"""<img src={gravatar_url} height={size} width={size}/>""")


@register.inclusion_tag('widgets/tags_banner.html', takes_context=True)
def tags_banner(context, limit=5, listing=False):

    request = context["request"]
    page = request.GET.get("page")

    tags = Post.objects.order_by("-pk").values("tags__name").annotate(Count('tags__name'))

    if listing:
        # Get the page info
        paginator = Paginator(tags, settings.TAGS_PER_PAGE)
        all_tags = paginator.get_page(page)
    else:
        all_tags = tags

    return dict(tags=all_tags, limit=limit, listing=listing, request=request)


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree, form):

    "Renders the post body"
    request = context['request']
    return dict(post=post, user=user, tree=tree, request=request,
                form=form,
                redir_field_name=const.REDIRECT_FIELD_NAME)


@register.inclusion_tag('widgets/subs_actions.html')
def subs_actions(post, user):

    if user.is_anonymous:
        sub = None
    else:
        sub = post.subs.filter(user=user).first()

    sub_type = models.Subscription.NO_MESSAGES if not sub else sub.type

    initial = dict(subtype=sub_type)

    form = forms.SubsForm(user=user, post=post, initial=initial)
    unsubbed = sub_type == models.Subscription.NO_MESSAGES

    button = "Follow" if unsubbed else "Update"

    return dict(post=post, form=form, button=button)


@register.inclusion_tag("widgets/forum_top_actionbar.html", takes_context=True)
def forum_top_actionbar(context, base_url="", objs=None):

    bar_objs = context.get("objs", objs)
    extra_context = dict(base_url=reverse(base_url), objs=bar_objs)
    context.update(extra_context)
    return context


@register.filter
def get_last_login(user):
    if user.profile.last_login:
        return f"visited {time_ago(user.profile.last_login)}"
    return f"visited {time_ago(user.profile.date_joined)}"


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

    color = "modcolor" if user.profile.is_moderator else ""

    if user.profile.score > 150:
        icon = f'<i class="ui bolt icon {color}"></i>'
    else:
        icon = f'<i class="ui genderless icon {color}"></i>'

    return mark_safe(icon)


@register.filter
def show_score(score):

    score = (score * 14) + 1
    return score


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


@register.filter
def time_ago(date):
    return util.time_ago(date=date)


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
        style = "gold"
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


@register.simple_tag
def render_comments(request, tree, post, next_url, project_uid=None,
                    comment_template='widgets/comment_body.html'):

    if post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree,
                                 comment_template=comment_template,
                                 next_url=next_url, project_uid=project_uid)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, comment_template, next_url,
                      project_uid=None):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(comment_template)
    comment_url = reverse("post_comment")

    def traverse(node):
        vote_url = reverse("vote")

        data = ['<div class="ui comment segments">']
        cont = {"post": node, 'user': request.user, 'request': request, "comment_url":comment_url,
                "vote_url": vote_url, "next_url":next_url, "redir_field_name":const.REDIRECT_FIELD_NAME,
                "project_uid": project_uid}
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


