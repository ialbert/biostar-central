import logging
import hashlib
import urllib.parse
from datetime import timedelta, datetime
from django.utils.timezone import utc
from django import template
from django.utils.safestring import mark_safe
from django.core.paginator import Paginator
from django.conf import settings
from django.contrib.auth import get_user_model

from biostar.shortcuts import reverse
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
def pages(objs, request=None, topic="", url="post_list_topic", uid=None):

    topic = request.GET.get('topic', topic)

    if topic:
        url = reverse(url, request=request, kwargs=dict(topic=topic))
    elif uid:
        url = reverse(url, request=request, kwargs=dict(uid=uid))
    else:
        url = ''

    return dict(objs=objs, url=url)


@register.inclusion_tag("widgets/message_menu.html")
def message_menu(inbox=None, unread=None, mentioned=None,
                 projects=None, outbox=None, request=None):

    return dict(inbox=inbox,unread=unread,mentioned=mentioned,
                projects=projects, outbox=outbox,request=request,
                active_tab=const.ACTIVE_TAB, const_in=const.INBOX,
                const_out=const.OUTBOX, const_unread=const.UNREAD)


@register.simple_tag
def only_enable_forum():

    return settings.ONLY_ENABLE_FORUM


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


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree, form, include_userbox=True, comment_view="post_comment",
              sub_redir="post_view", vote_view="update_vote", sub_view="subs_action"):
    #TODO: this is really temporary ( the view names cannot be a string here.)
    "Renders the post body"
    request = context['request']

    vote_url = reverse(vote_view, request=request, kwargs=dict(uid=post.uid))
    sub_url = reverse(sub_view, request=request, kwargs=dict(uid=post.uid))
    next_url = reverse(sub_redir, request=request, kwargs=dict(uid=post.uid))
    comment_url = reverse(comment_view, request=request, kwargs=dict(uid=post.uid))

    return dict(post=post, user=user, tree=tree, request=request,
                form=form, include_userbox=include_userbox, comment_url=comment_url,
                vote_redir=sub_redir, sub_url=sub_url, vote_url=vote_url,
                comment_view=comment_view, next_url=next_url, vote_view=vote_view,
                redir_field_name=const.REDIRECT_FIELD_NAME)


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
def feed(user, post=None, limit=7):
    #TODO: temporary feed

    # Show similar posts when inside of a view
    if post:
        return
    post_set = Post.objects.exclude(status=Post.DELETED).all()

    recent_votes = Vote.objects.filter(type=Vote.UP)[:limit]
    # Needs to be put in context of posts
    recent_votes = post_set.filter(votes__in=recent_votes).order_by("?")

    # TODO:change
    recent_locations = User.objects.filter(post__in=post_set).order_by("?").distinct()

    recent_locations = set([x for x in recent_locations if x.profile.location][:limit])

    recent_awards = ''
    recent_replies = post_set.filter(type__in=[Post.COMMENT, Post.ANSWER],
                                     ).order_by("?")[:limit]

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
def listing(posts=None, messages=None):

    is_post = True if posts else False
    is_messages = True if messages else False

    objs = posts or messages

    return dict(is_post=is_post, is_messages=is_messages, objs=objs)



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
    if post.has_accepted:
        style = "accepted"
    elif post.reply_count > 0:
        style = "answered"
    elif post.comment_count > 0:
        style = "commented"
    else:
        style = "unanswered"

    return style


@register.simple_tag
def render_comments(request, post, comment_view, vote_view, next_url, comment_template='widgets/comment_body.html'):

    user = request.user
    thread = Post.objects.get_thread(post.parent, user)
    # Build tree
    tree = auth.build_tree(thread=thread, tree={})

    if tree and post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree,
                                 comment_template=comment_template, comment_view=comment_view,
                                 vote_view=vote_view, next_url=next_url)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, comment_template, comment_view, vote_view, next_url):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(comment_template)

    def traverse(node):
        comment_url = reverse(comment_view, request=request, kwargs=dict(uid=node.uid))
        vote_url = reverse(vote_view, request=request, kwargs=dict(uid=node.uid))

        data = ['<div class="comment">']
        cont = {"post": node, 'user': request.user, 'request': request, "comment_url":comment_url,
                "vote_url":vote_url, "next_url":next_url, "redir_field_name":const.REDIRECT_FIELD_NAME}
        html = body.render(cont)
        data.append(html)
        for child in tree.get(node.id, []):

            data.append('<div class="comments" >')
            data.append(traverse(child))
            data.append("</div>")

        data.append("</div>")

        return '\n'.join(data)

    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:

        coll.append(traverse(node))
    return '\n'.join(coll)


