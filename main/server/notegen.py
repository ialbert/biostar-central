"""
This module can generate messages
"""
from main.server.const import *

def userlink(user, short=False):
    if short:
        return '%s' % user.profile.display_name
    else:
        return '[%s](%s)' % (user.profile.display_name, user.profile.get_absolute_url() )

def postlink(post, size=35):
    root  = post.root or post
    title = root.title
    title = title if len(title) < size else '%s...' %  title[:size]
    return '[%s: %s]( %s)' % (post.get_type_display(), title, post.get_absolute_url())

def badgelink(badge):
    return '[%s](%s)' % (badge.name, badge.get_absolute_url())

def post_moderator_action(user, post):
    action = post.get_status_display()
    text   = '%s set status to %s on %s' % (userlink(user), action, postlink(post))
    return text

def user_moderator_action(user, target):
    action = target.profile.get_status_display()
    text   = '%s set status %s on %s' % (userlink(user), action, userlink(target))
    return text

def chop(text, size):
    #text = text.encode('ascii', errors='replace')
    return text if len(text) < size else '%s...' % text[:size]
    
def post_action(user, post, size=250):
    text   = '%s by %s : %s' % (postlink(post), userlink(user, short=False), chop(post.content, size=size))
    return text

def awardnote(badge):
    text = "Congratulations! You've just earned the **%s** badge!" % badgelink(badge)
    return text

def private_message(user, target, text):
    text = "%s by %s" % (userlink(user), text)
    return text