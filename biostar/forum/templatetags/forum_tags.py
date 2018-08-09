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

from biostar.accounts.models import Profile
from biostar.utils.shortcuts import reverse
from biostar.forum.models import Post, Vote, Message
from biostar.forum import auth, forms, models, const


User = get_user_model()

logger = logging.getLogger("engine")

register = template.Library()

def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.inclusion_tag('widgets/user_box.html')
def user_box(user):

    return dict(user=user)


@register.inclusion_tag('widgets/pages.html')
def pages(objs, request):

    topic = request.GET.get('topic', request.GET.get("active"))

    url = request.path

    return dict(objs=objs, url=url, topic=topic)


@register.inclusion_tag("widgets/message_menu.html")
def message_menu(extra_tab=None, request=None):

    extra = {extra_tab: "active"}
    context = dict(request=request, active_tab=const.ACTIVE_MESSAGE_TAB,
                   const_in=const.INBOX, const_out=const.OUTBOX,
                   const_unread=const.UNREAD)
    context.update(extra)
    return context


@register.inclusion_tag('widgets/forum_menubar.html', takes_context=True)
def forum_menubar(context, request=None):
    user = context.request.user

    return dict(user=user, request=request)


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
        'd': 'identicon',
    }
    )

    return mark_safe(f"""<img src={gravatar_url} height={size} width={size}/>""")


@register.inclusion_tag('widgets/tags_banner.html', takes_context=True)
def tags_banner(context, limit=7, listing=False):

    request = context["request"]
    page = request.GET.get("page")
    default = ["job", "news"]

    default = list(map(lambda x: dict(tags__name=x, tags__name__count=1), default))

    tags = list(Post.objects.values("tags__name").annotate(Count('tags__name')))[:limit]

    tags = list(filter(lambda x: x["tags__name"] is not None, tags))

    if listing:
        # Get the page info
        paginator = Paginator(tags, settings.TAGS_PER_PAGE)
        all_tags = paginator.get_page(page)
    else:
        all_tags = default + tags

    return dict(tags=all_tags, limit=limit, listing=listing, request=request)


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree, form, include_userbox=True, comment_url="/",
              sub_redir="post_view", vote_view="update_vote", sub_view="subs_action", project_uid=None):

    "Renders the post body"
    request = context['request']

    vote_url = reverse(vote_view, request=request, kwargs=dict(uid=post.uid))
    sub_url = reverse(sub_view, request=request, kwargs=dict(uid=post.uid))
    next_url = reverse(sub_redir, request=request, kwargs=dict(uid=post.uid))

    return dict(post=post, user=user, tree=tree, request=request,
                form=form, include_userbox=include_userbox, comment_url=comment_url,
                vote_redir=sub_redir, sub_url=sub_url, vote_url=vote_url, next_url=next_url, vote_view=vote_view,
                redir_field_name=const.REDIRECT_FIELD_NAME, project_uid=project_uid)


@register.inclusion_tag('widgets/subs_actions.html')
def subs_actions(post, user, next, sub_url):

    if user.is_anonymous:
        sub = None
    else:
        sub = models.Subscription.get_sub(user=user, post=post).first()

    sub_type = models.Subscription.NO_MESSAGES if not sub else sub.type

    initial = dict(subtype=sub_type)

    form = forms.SubsForm(user=user, post=post, initial=initial)
    unsubbed = sub_type == models.Subscription.NO_MESSAGES

    button = "Follow" if unsubbed else "Update"

    return dict(post=post, form=form, button=button, next=next, sub_url=sub_url,
                redir_field_name=const.REDIRECT_FIELD_NAME)


@register.filter
def show_email(user):

    try:
        head, tail = user.email.split("@")
        email = head[0] + "*" * 10 + tail
    except:
        return user.email[0] + "*" * 10

    return email


@register.inclusion_tag('widgets/feed.html')
def feed(user, post=None):
    #TODO: temporary feed

    # Show similar posts when inside of a view
    if post:
        return

    recent_votes = Vote.objects.filter(type=Vote.UP)[:settings.VOTE_FEED_COUNT]
    # Needs to be put in context of posts
    recent_votes = recent_votes.select_related("post")

    recent_locations = User.objects.filter(
        ~Q(profile__location="")).select_related("profile").distinct()[:settings.LOCATION_FEED_COUNT]

    recent_awards = ''
    recent_replies = Post.objects.filter(~Q(status=Post.DELETED), type__in=[Post.COMMENT, Post.ANSWER]
                                     ).select_related("author__profile", "author")[:settings.REPLIES_FEED_COUNT]

    context = dict(recent_votes=recent_votes, recent_awards=recent_awards,
                   recent_locations=recent_locations, recent_replies=recent_replies,
                   post=post, user=user)

    return context


@register.filter
def show_score_icon(score):

    icon = "small circle"
    if score > 500:
        icon = "small star"

    score_icon = f'<i class="ui {icon} icon"></i>'

    return mark_safe(score_icon)


@register.filter
def show_score(score):

    score = (score * 14) + 1
    return score


@register.inclusion_tag('widgets/user_info.html')
def user_info(post, by_diff=False):

    return dict(post=post, by_diff=by_diff)


@register.simple_tag
def get_posts(user, request, per_page=20):

    posts = Post.objects.my_posts(target=user, user=user)
    page = request.GET.get("page", 1)

    paginator = Paginator(posts, per_page=per_page)
    page = page if page is not None else 1

    objs = paginator.get_page(page)

    return objs


@register.inclusion_tag('widgets/listing.html')
def listing(posts=None, messages=None, discussion_view=False):

    is_post = True if posts else False
    is_messages = True if messages else False

    objs = posts or messages

    return dict(is_post=is_post, is_messages=is_messages, objs=objs, discussion_view=discussion_view)


@register.simple_tag
def get_top_padding(post):
    #TODO: temporary solve

    if len(post.get_title()) >= 63:
        return "small-padding"
    return ""


@register.filter
def show_nonzero(value):
    "The purpose of this is to return value or empty"
    return value if value else ''


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


@register.simple_tag
def object_count(request, otype):

    user = request.user

    if user.is_authenticated:
        if otype == "post":
            return Post.objects.my_posts(target=user, user=user).count()
        if otype == "follow":
            # Stuff that produces notifications
            query = models.Subscription.objects.exclude(type=models.Subscription.NO_MESSAGES).filter(user=user)
            return query.count()
        if otype == "bookmark":
            return Post.objects.my_bookmarks(user).count()
        if otype == "votes":
            return  Post.objects.my_post_votes(user).distinct().count()
        if otype == "message":
            # Return the count stored in the message
            return user.profile.new_messages
        if otype == "unread":
            return Message.objects.filter(recipient=user, unread=True).count()
        if otype =="inbox":
            return Message.objects.inbox_for(user=user).count()
        if otype == "outbox":
            return Message.objects.outbox_for(user=user).count()

    return 0



@register.filter
def preview_message(text, limit=130):

    return text if len(text) <= limit else text[:limit] + " ..."



@register.simple_tag
def vote_icon(user, post, vtype):

    main_map = {"bookmark":{"icon":"bookmark", "vote":Vote.BOOKMARK},
                "upvote":{"icon":"thumbs up", "vote":Vote.UP}
                }

    icon, vote_type = main_map[vtype]["icon"], main_map[vtype]["vote"]

    if user.is_authenticated:
        vote = Vote.objects.filter(author=user, post=post, type=vote_type).first()
    else:
        vote = None

    msg = f"{icon} icon"

    if not vote:
        msg += " outline"

    return mark_safe(msg)


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

    if post.type == Post.JOB:
        style = "orange"
    elif post.type == Post.TUTORIAL:
        style = "blue"
    elif post.type == Post.TOOL:
        style = "darkgreen"
    elif post.type == Post.FORUM:
        style = "gold"
    elif post.type == Post.NEWS:
        style = "purple"
    elif post.has_accepted:
        style = "olive"
    elif post.reply_count > 0:
        style = "teal"
    elif post.comment_count > 0:
        style = "grey"
    else:
        style = "maroon"

    return style


@register.simple_tag
def render_comments(request, post, comment_url, vote_view, next_url, project_uid=None,
                    comment_template='widgets/comment_body.html'):

    user = request.user
    thread = Post.objects.get_thread(post.parent, user)
    # Build tree
    tree = auth.build_tree(thread=thread, tree={})

    if tree and post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree,
                                 comment_template=comment_template, comment_url=comment_url,
                                 vote_view=vote_view, next_url=next_url, project_uid=project_uid)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, comment_url, comment_template, vote_view, next_url,
                      project_uid=None):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(comment_template)

    def traverse(node):
        vote_url = reverse(vote_view, request=request, kwargs=dict(uid=node.uid))

        data = ['<div class="ui segments">']
        cont = {"post": node, 'user': request.user, 'request': request, "comment_url":comment_url,
                "vote_url":vote_url, "next_url":next_url, "redir_field_name":const.REDIRECT_FIELD_NAME,
                "project_uid": project_uid}
        html = body.render(cont)
        data.append(html)
        for child in tree.get(node.id, []):

            data.append('<div class="ui segment comment" >')
            data.append(traverse(child))
            data.append("</div>")

        data.append("</div>")

        return '\n'.join(data)

    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:

        coll.append(traverse(node))
    return '\n'.join(coll)


