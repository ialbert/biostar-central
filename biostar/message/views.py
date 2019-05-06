from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.conf import settings
from django.core.paginator import Paginator
from biostar.accounts.models import Profile

from .models import Message
from .decorators import message_access
from .const import *
from . import auth, forms


def get_messages(user, listing):

    message_map = dict(
        inbox=Message.objects.filter(recipient=user),
        unread=Message.objects.filter(recipient=user, unread=True),
        outbox=Message.objects.filter(sender=user),
        mentioned=Message.objects.filter(recipient=user, source=Message.MENTIONED)
    )

    # Get the messages based on queried topic
    messages = message_map.get(listing, message_map["inbox"])

    messages = messages.select_related("recipient", "sender", "sender__profile",
                                       "recipient__profile")
    messages = messages.order_by("-sent_date")
    return messages


def message_list(request, template="message_list.html", listing=INBOX):

    page = request.GET.get("page", 1)

    messages = get_messages(user=request.user, listing=listing)

    # Get the pagination info
    paginator = Paginator(messages, settings.MESSAGES_PER_PAGE)
    messages = paginator.get_page(page)

    context = dict(message="active", all_messages=messages, extra_tab_name=listing)

    user = request.user

    Profile.objects.filter(user=user).update(new_messages=0)

    return render(request, template, context)


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


@message_access(access_to=INBOX)
@login_required
def inbox_view(request, uid):
    "Checks to see you are the recipient to a message."

    return message_view(request, uid)


@message_access(access_to=OUTBOX)
@login_required
def outbox_view(request, uid):
    "Checks to see if you are the sender of the message."

    return message_view(request, uid)


@login_required
def inbox_list(request):

    # Get the unread messages from inbox
    listing = request.GET.get("active", INBOX)

    return message_list(request, template="message_list.html", listing=listing)


@login_required
def outbox_list(request):
    return message_list(request, template="message_list.html", listing=OUTBOX)


@login_required
def message_compose(request):

    # Get the message author
    author = request.user

    # Load the form.

    form = forms.Compose()

    if request.method == "POST":
        form = forms.Compose(data=request.POST)
        if form.is_valid():
            form.save()
            1/0
            pass

    context = dict(form=form)

    return render(request, "message_compose.html", context)

