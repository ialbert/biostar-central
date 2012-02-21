"""
This module can generate messages
"""
from main.server.const import *

def userlink(user):
    return '[%s](%s)' % (user.profile.display_name, user.profile.get_absolute_url() )

def postlink(post, size=35):
    root  = post.root or post
    title = root.title if len(root.title) < size else '%s...' %  root.title[:size] 
    return '[%s](%s)' % (title, root.get_absolute_url())

def badgelink(badge):
    return '[%s](%s)' % (badge.name, badge.get_absolute_url() )

def post_moderator_action(user, post, action):
    action = REV_ACTION_MAP.get(action, '???')
    text  = '%s used %s on %s' % (userlink(user), action, postlink(post))
    return text

def post_action(user, post):
   
    post_type = int(post.type)

    if post_type == POST_QUESTION:
        action = 'asked'
    elif post_type == POST_ANSWER:
        action = 'answered'
    elif post_type == POST_COMMENT:
        action = 'commented on'
    else:
        action = "posted"
    
    size = 250
    if len(post.content)<size:
        content = post.content
    else:
        content = '%s...' % post.content[:size]
    
    text   = '%s %s %s %s' % (userlink(user), action, postlink(post), content)
    return text

def suspend(user):
    text = "suspended %s" % userlink(user)
    return text

def reinstate(user):
    text = "reinstated %s" % userlink(user)
    return text

def badgenote(badge):
    text = "Congratulations! You've just earned the **%s** badge!" % badgelink(badge)
    return text