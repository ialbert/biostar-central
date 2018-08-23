from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.conf import settings
from django.core.paginator import Paginator
from biostar.accounts.models import Profile

from .models import Message
from .decorators import message_access
from .const import *
from . import auth


def message_view(request, uid):

    base_message = Message.objects.filter(uid=uid).first()

    # Build the message tree from bottom up
    tree = auth.build_msg_tree(msg=base_message, tree=[])

    # Update the unread flag
    Message.objects.filter(pk=base_message.pk).update(unread=False)

    active_tab = request.GET.get(ACTIVE_MESSAGE_TAB, INBOX)

    context = dict(base_message=base_message, tree=tree, extra_tab="active",
                   extra_tab_name=active_tab)

    return render(request, "message_view.html", context=context)


def message_list(request, template="message_list.html", active_tab=None, listing=INBOX):

    page = request.GET.get("page", 1)

    objs = auth.list_message_by_topic(request=request, topic=active_tab).order_by("-pk")

    # Get the page info
    paginator = Paginator(objs, settings.MESSAGES_PER_PAGE)
    objs = paginator.get_page(page)

    context = dict(is_inbox=listing == INBOX, field_name=ACTIVE_MESSAGE_TAB,
                   extra_tab_name=active_tab.capitalize(),extra_tab="active", objs=objs)

    context.update({active_tab: ACTIVE_MESSAGE_TAB})

    user = request.user

    Profile.objects.filter(user=user).update(new_messages=0)

    return render(request, template, context)


@message_access(access_to=INBOX)
@login_required
def inbox_message_view(request, uid):
    "Checks to see you are the recipient to a message."

    return message_view(request, uid)


@message_access(access_to=OUTBOX)
@login_required
def outbox_message_view(request, uid):
    "Checks to see if you are the sender of the message."

    return message_view(request, uid)


@login_required
def inbox_list(request):

    active_tab = request.GET.get(ACTIVE_MESSAGE_TAB, INBOX)
    active_tab = active_tab if (active_tab in MESSAGE_TABS) else INBOX

    return message_list(request, template="message_list.html", active_tab=active_tab,
                        listing=INBOX)


@login_required
def outbox_list(request):
    return message_list(request, template="message_list.html", active_tab=OUTBOX,
                        listing=OUTBOX)


@login_required
def message_compose(request):

    return

