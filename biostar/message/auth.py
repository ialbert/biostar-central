from .models import Message
from biostar.forum.util import fixcase
from . import const


def query_topic(user, topic):

    mapper = {
        const.INBOX: dict(func=Message.objects.inbox_for, params=dict(user=user)),
        const.UNREAD: dict(func=Message.objects.filter, params=dict(recipient=user, unread=True)),
        const.OUTBOX: dict(func=Message.objects.outbox_for, params=dict(user=user)),
        const.MENTIONED: dict(func=Message.objects.inbox_for, params=dict(user=user),
                          apply=lambda q: q.filter(source=Message.MENTIONED)),
        }

    if mapper.get(topic):
        func, params = mapper[topic]["func"], mapper[topic].get("params")
        apply_extra = mapper[topic].get("apply", lambda q: q)
        query = apply_extra(func(**params))
    else:
        query = None

    return query


def list_message_by_topic(request, topic):

    user = request.user
    # One letter tags are always uppercase
    topic = fixcase(topic)

    query = query_topic(user=user, topic=topic)

    if query is None:
        query = Message.objects.inbox_for(user=user)

    return query


def build_msg_tree(msg, tree=[]):
    "Build a flat tree with the msg being at the base ( index=0)"

    # Add the current message
    tree.append(msg)

    # Check if it has a parent
    # and recursively add that to that tree.
    if msg.parent_msg and msg.parent_msg != msg:
        tree.append(msg.parent_msg)
        build_msg_tree(msg=msg.parent_msg, tree=tree)

    # End of tree, at the root message
    return tree


def create_messages(body, sender, recipient_list, subject="", parent=None,
                    source=Message.REGULAR, mtype=None):
    "Create batch message from sender for a given recipient_list"

    subject = subject or f"Message from : {sender.profile.name}"

    msg_list = []

    #TODO: do a bulk create for the whole recipeint list.
    for rec in recipient_list:

        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject,
                                     body=body, parent_msg=parent, type=mtype, source=source)
        msg_list.append(msg)

    return msg_list