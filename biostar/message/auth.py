import uuid
from datetime import datetime

from django.utils.timezone import utc

from .models import Message
from biostar.forum.util import fixcase
from biostar.accounts.models import Profile
from . import const

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

def now():
    return datetime.utcnow().replace(tzinfo=utc)


def build_msg_tree(msg, tree=[]):
    "Build a flat tree with the msg being at the base ( index=0)"

    # Add the current message
    tree.append(msg)

    # Check if it has a parent
    # and recursively add that to the tree.
    if msg.parent_msg and msg.parent_msg != msg:
        tree.append(msg.parent_msg)
        build_msg_tree(msg=msg.parent_msg, tree=tree)

    # End of tree, at the root message
    return tree


def create_local_messages(body, sender, rec_list, subject="", parent=None,
                          source=Message.REGULAR, mtype=Profile.LOCAL_MESSAGE):
    "Create batch message from sender for a given recipient_list"

    subject = subject or f"Message from : {sender.profile.name}"
    msgs = []
    for rec in rec_list:
        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject, source=source,
                                     sent_date=now(), uid=get_uuid(10), body=body, parent_msg=parent, type=mtype)

        msgs.append(msg)

    return msgs