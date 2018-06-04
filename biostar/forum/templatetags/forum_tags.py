import logging
from datetime import timedelta, datetime
from django.utils.timezone import utc
from django import template
from django.utils.safestring import mark_safe

from biostar.forum.models import Post
from biostar.forum import auth, forms, models


logger = logging.getLogger("engine")

register = template.Library()

def now():
    return datetime.utcnow().replace(tzinfo=utc)


@register.inclusion_tag('widgets/post_body.html', takes_context=True)
def post_body(context, post, user, tree, form, add_comment=False):
    "Renders the post body"

    return dict(post=post, user=user, tree=tree, request=context['request'],
                add_comment=add_comment, form=form)


@register.inclusion_tag('widgets/subs_actions.html')
def subs_actions(post, user):


    if user.is_anonymous:
        sub = None
    else:
        sub = models.Subscription.get_sub(user=user, post=post).first()

    sub_type = models.Subscription.NO_MESSAGES if not sub else sub.type

    initial = dict(subtype=sub_type)

    form = forms.SubsForm(user=user, post=post, initial=initial)
    unsubbed = sub_type == models.Subscription.NO_MESSAGES

    button = "Follow" if unsubbed else "Update"

    return dict(post=post, form=form, button=button)


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

    return 0



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
def render_comments(request, post, comment_template='widgets/comment_body.html'):


    user = request.user
    thread = Post.objects.get_thread(post.parent, user)
    # Build tree
    tree = auth.build_tree(thread=thread, tree={})

    if tree and post.id in tree:
        text = traverse_comments(request=request, post=post, tree=tree,
                                 comment_template=comment_template)
    else:
        text = ''

    return mark_safe(text)


def traverse_comments(request, post, tree, comment_template):
    "Traverses the tree and generates the page"

    body = template.loader.get_template(comment_template)

    def traverse(node):
        data = ['<div class="comment">']
        cont = {"post": node, 'user': request.user, 'request': request}
        html = body.render(cont)
        data.append(html)
        for child in tree.get(node.id, []):

            data.append('<div class="comments">')
            data.append(traverse(child))
            data.append("</div>")

        data.append("</div>")

        return '\n'.join(data)

    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:

        coll.append(traverse(node))
    return '\n'.join(coll)


