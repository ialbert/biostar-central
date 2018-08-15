import logging
from django import template
from biostar.message import const, auth, models


logger = logging.getLogger("engine")

register = template.Library()


@register.simple_tag
def get_all_message_count(request):

    user = request.user
    outbox = inbox = projects = mentioned = unread = 0

    if user.is_authenticated:
        inbox = models.Message.objects.inbox_for(user=user)
        outbox = models.Message.objects.outbox_for(user=user).count()
        unread = inbox.filter(unread=True).count()

    context = dict(outbox=outbox, inbox=inbox.count(), projects=projects, mentioned=mentioned,
                   unread=unread)

    return context


@register.simple_tag
def message_count(request, otype):

    user = request.user
    count = 0

    if user.is_authenticated:

        if otype == "message":
            count = user.profile.new_messages
        else:
            query = auth.query_topic(user=user, topic=otype)
            count = count if query is None else query.count()

    return count


@register.inclusion_tag("widgets/message_menu.html")
def message_menu(extra_tab=None, request=None):

    extra = {extra_tab: "active"}
    context = dict(request=request, active_tab=const.ACTIVE_MESSAGE_TAB,
                   const_in=const.INBOX, const_out=const.OUTBOX,
                   const_unread=const.UNREAD)
    context.update(extra)
    return context
