import logging
from django import template
from biostar.message import const, models

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
def get_all_message_count(request):

    user = request.user
    outbox = inbox = projects = mentioned = unread = 0

    if user.is_authenticated:
        inbox = models.Message.objects.inbox_for(user=user)
        outbox = models.Message.objects.outbox_for(user=user).count()
        unread = inbox.filter(unread=True).count()
        mentioned = inbox.filter(source=models.Message.MENTIONED).count()

    context = dict(outbox=outbox, inbox=inbox.count(), projects=projects, mentioned=mentioned,
                   unread=unread)

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


@register.inclusion_tag("widgets/message_top_actionbar.html", takes_context=True)
def message_top_actionbar(context, with_pages=False):
    extra_context = dict(with_pages=with_pages)
    context.update(extra_context)
    return context


@register.simple_tag
def is_unread(user, message):

    if message.recipient == user and message.unread:
        return "unread-message"

    return ""


@register.simple_tag
def is_inbox(message, target):

    return message.recipient == target


@register.inclusion_tag("widgets/message_menu.html")
def message_menu(tab=None, request=None):

    extra = {tab: "active"}
    context = dict(request=request, active_tab=const.ACTIVE_MESSAGE_TAB,
                   const_in=const.INBOX, const_out=const.OUTBOX,
                   const_unread=const.UNREAD)
    context.update(extra)
    return context
