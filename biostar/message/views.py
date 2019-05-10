from datetime import timedelta
from django.shortcuts import render, redirect, reverse
from django.contrib.auth.decorators import login_required
from django.conf import settings
from django.contrib import messages
from django.core.paginator import Paginator
from biostar.accounts.models import Profile

from .models import Message
from .decorators import message_access
from .const import *
from . import auth, forms


LIMIT_MAP = dict(
    all=0,
    today=1,
    week=7,
    month=30,
    year=365
)

ORDER_MAPPER = dict(
    rsent="sent_date",
    sent="-sent_date",
    rep="sender__profile__score",

)


def get_messages(user, listing, limit=0, order="sent"):

    msg_map = dict(
        inbox=Message.objects.filter(recipient=user),
        unread=Message.objects.filter(recipient=user, unread=True),
        outbox=Message.objects.filter(sender=user),
        mentioned=Message.objects.filter(recipient=user, source=Message.MENTIONED)
    )

    # Get the messages based on queried topic
    msgs = msg_map.get(listing, msg_map["inbox"])

    if ORDER_MAPPER.get(order):
        ordering = ORDER_MAPPER.get(order)
        msgs = msgs.order_by(ordering)
    else:
        msgs = msgs.order_by("-sent_date")

    days = LIMIT_MAP.get(limit, 0)
    # Apply time limit if required.
    if days:
        delta = auth.now() - timedelta(days=days)
        msgs = msgs.filter(sent_date__gt=delta)

    msgs = msgs.select_related("recipient", "sender", "sender__profile",
                               "recipient__profile")
    #msgs = msgs.order_by("-sent_date")
    return msgs


def message_list(request, template="message_list.html", listing=INBOX):

    page = request.GET.get("page", 1)
    limit = request.GET.get("limit", 0)
    order = request.GET.get("order", "sent")
    messages = get_messages(user=request.user, listing=listing, limit=limit, order=order)

    # Get the pagination info
    paginator = Paginator(messages, settings.MESSAGES_PER_PAGE)
    messages = paginator.get_page(page)

    context = dict(message="active", all_messages=messages, order=order, limit=limit,
                   extra_tab_name=listing)

    user = request.user

    Profile.objects.filter(user=user).update(new_messages=0)

    return render(request, template, context)


def message_view(request, msg):

    # Build the message tree from bottom up
    tree = auth.build_msg_tree(msg=msg, tree=[])

    active_tab = request.GET.get(ACTIVE_MESSAGE_TAB, INBOX)

    context = dict(base_message=msg, tree=tree, extra_tab="active",
                   extra_tab_name=active_tab)

    return render(request, "message_view.html", context=context)


@message_access(access_to=INBOX)
@login_required
def inbox_view(request, uid):
    "Checks to see you are the recipient to a message."

    msg = Message.objects.filter(uid=uid).first()
    # Update the unread flag
    Message.objects.filter(pk=msg.pk).update(unread=False)

    return message_view(request, msg=msg)


@message_access(access_to=OUTBOX)
@login_required
def outbox_view(request, uid):
    "Checks to see if you are the sender of the message."
    msg = Message.objects.filter(uid=uid).first()
    return message_view(request, msg=msg)


@login_required
def inbox_list(request):

    # Get the unread messages from inbox
    listing = request.GET.get("active", INBOX)

    return message_list(request, template="message_list.html", listing=listing)


@login_required
def outbox_list(request):
    return message_list(request, template="message_list.html", listing=OUTBOX)


def block_user(request):
    return


def report_spam(request):
    return


@login_required
def message_compose(request):

    # Get the message author
    author = request.user

    # Get an initial user to add to recipients list
    inital_user = request.GET.get("initial", "")
    initial = dict(recipients=inital_user)
    # Load the form.
    form = forms.Compose(initial=initial)

    if request.method == "POST":
        form = forms.Compose(data=request.POST)
        if form.is_valid():
            form.save(sender=author)
            messages.success(request, "Sent message to recipients")
            return redirect(reverse("outbox"))

    context = dict(form=form, compose='active')

    return render(request, "message_compose.html", context)

