from datetime import timedelta
from django.shortcuts import render, redirect, reverse
from django.contrib.auth.decorators import login_required
from django.conf import settings
from django.contrib import messages
from django.core.paginator import Paginator
from biostar.accounts.models import Profile
from biostar.utils.decorators import object_exists

from .models import Message
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


def index(request):
    context = dict()
    return render(request, "message/index.html", context)


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
    return msgs

@login_required
def message_list(request, listing="inbox"):

    user = request.user

    page = request.GET.get("page", 1)
    limit = request.GET.get("limit", 0)
    order = request.GET.get("order", "sent")

    listing = request.GET.get("active", None) or listing
    #print(request.GET.get("active", None))

    messages = get_messages(user=request.user, listing=listing, limit=limit, order=order)

    # Get the pagination info
    paginator = Paginator(messages, settings.MESSAGES_PER_PAGE)
    messages = paginator.get_page(page)

    context = dict(tab="messages", all_messages=messages, order=order, limit=limit,
                   tab_name=listing)

    Profile.objects.filter(user=user).update(new_messages=0)

    return render(request, "message/message_list.html", context)


@login_required
def outbox(request):
    return message_list(request, listing="outbox")


@login_required
def reply(request, uid):
    parent_msg = Message.objects.filter(uid=uid).first()
    form = forms.Reply()
    context = dict(msg=parent_msg, form=form)
    return render(request, "message/reply.html", context=context)


@object_exists(klass=Message)
@login_required
def inbox_view(request, uid):
    "Checks to see you are the recipient to a message."

    msg = Message.objects.filter(recipient=request.user, uid=uid).first()
    # Update the unread flag
    Message.objects.filter(pk=msg.pk).update(unread=False)

    # Get the reply form
    if request.method == "POST":
        form = forms.Reply(data=request.POST)
        if form.is_valid():
            uid = form.save(sender=request.user, recipients=msg.sender.username, parent=msg)
            # Redirect to outbox of recently sent with message
            return redirect(reverse("outbox_view", kwargs=dict(uid=uid)))

    # Build the message tree from bottom up
    tree = auth.build_msg_tree(msg=msg, tree=[])

    active_tab = request.GET.get("active", INBOX)

    context = dict(base_message=msg, tree=tree, tab_name=active_tab, message=msg, tab="messages")
    return render(request, "message/message_view.html", context=context)


@object_exists(klass=Message)
@login_required
def outbox_view(request, uid):
    "Checks to see if you are the sender of the message."
    sender = request.user
    msg = Message.objects.filter(sender=sender, uid=uid).first()
    if not msg:
        messages.error(request, "Message does not exist in your outbox")
        return redirect(reverse("outbox"))

    # Build the message tree from bottom up
    tree = auth.build_msg_tree(msg=msg, tree=[])

    context = dict(base_message=msg, tree=tree, tab_name=OUTBOX, message=msg, tab="messages")
    return render(request, "message/message_view.html", context=context)


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

    context = dict(form=form, tab_name="compose", tab="messages")

    return render(request, "message/message_compose.html", context)

