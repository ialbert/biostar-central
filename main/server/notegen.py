"""
This module can generate messages
"""
from main.server.const import *

def moderator_action(user, post, action):
    action = REV_ACTION_MAP.get(action, 'undefined action?')
    name  = '<a href="%s">%s</a>' % (user.profile.get_absolute_url(), user.profile.display_name)
    title = '<a href="%s">%s</a>' % (post.get_absolute_url(), post.title)
    text  = '%s used %s on %s' % (name, action, title)
    return text
