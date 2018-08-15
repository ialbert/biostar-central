from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.conf import settings
from django.core.paginator import Paginator
from biostar.accounts.models import Profile
from biostar.forum.decorators import object_exists

from .models import Message
from .decorators import message_access
from .const import *
from . import auth


@object_exists(klass=Message, url="message_list")
@message_access
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


@login_required
def message_list(request):

    active_tab = request.GET.get(ACTIVE_MESSAGE_TAB, INBOX)

    active_tab = active_tab if (active_tab in MESSAGE_TABS) else INBOX

    page = request.GET.get("page", 1)

    objs = auth.list_message_by_topic(request=request, topic=active_tab).order_by("-pk")

    # Get the page info
    paginator = Paginator(objs, settings.MESSAGES_PER_PAGE )
    objs = paginator.get_page(page)

    context = {active_tab: ACTIVE_MESSAGE_TAB, "not_outbox": active_tab != OUTBOX,
               "field_name": ACTIVE_MESSAGE_TAB, "extra_tab_name": active_tab.capitalize(),
               "extra_tab": "active", "objs": objs}

    user = request.user

    Profile.objects.filter(user=user).update(new_messages=0)

    return render(request, "message_list.html", context)
