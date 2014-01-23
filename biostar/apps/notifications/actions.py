"""
Renders messages based on
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime

def new_post_created(sender, instance, created, *args, **kwargs):
    "When a new post is created send notifications to all subscribers"

    pass
    # Get all subscriptions for the root post.
    #for sub in instance.root.subs.exclude(user=instance.author):
    #    #print (sub.post.id, sub.user)
    #    pass