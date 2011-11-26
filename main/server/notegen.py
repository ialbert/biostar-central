"""
This module can generate messages
"""
from main.server.const import *

def userlink(user):
    return '[%s](%s)' % (user.profile.display_name, user.profile.get_absolute_url() )

def postlink(post):
    root = post.get_root()
    size = 30
    if len(root.title)<=size:
        title = root.title
    else:
        title = '%s...' %  root.title[:size] 
    return '[%s](%s%s/#%s/)' % (title, root.get_absolute_url(), root.slug, post.id)

def post_moderator_action(user, post, action):
    action = REV_ACTION_MAP.get(action, '???')
    text  = '%s used %s on %s' % (userlink(user), action, postlink(post))
    return text

def post_action(user, post):
    
    if post.post_type == POST_ANSWER:
        action = 'answered'
    elif post.post_type == POST_COMMENT:
        action = 'commented on'
    else:
        action = "did something to"
    size = 150
    if len(post.content)<size:
        content = post.content
    else:
        content = '%s...' % post.content[:size]
    text   = '%s %s %s with *%s*' % (userlink(user), action, postlink(post), content)
    return text