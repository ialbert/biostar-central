"""
User access authorization
""" 

import logging
logger = logging.getLogger(__name__)

def authorize_user_edit(target, user, strict=True):
    """
    Authorizes writing a target user by another user.  
    Also sets a new attribute 'writeable' on the post.
    Moderators may edit regular users, administrators may edit moderators and regular users.
    Strict mode raises an immediate exception.
    """
    if user.is_anonymous():
        writeable = False
    elif target == user:
        writeable = True
    elif target.profile.is_admin:
        # admins may not be moderated directly
        writeable = False  
    elif user.profile.is_admin:
        # admins may moderate everyone else
        writeable = True  
    elif target.profile.is_moderator:
        # at this point only moderators are left
        writeable = False        
    elif user.profile.is_moderator :
        # at this point the target is a regular user
        writeable = True
    else:
        # forbid access otherwise
        writeable = False
    
    if strict and not writeable:
        msg = 'user %s write access to user %s denied' % (user.id, target.id)
        logger.error(msg)
        raise Exception(msg)
    
    target.writeable = writeable
    return target.writeable

def authorize_post_edit(user, post, strict=True):
    """
    Authorizes editing a post by a user. Also sets a new attribute 'writeable' on the post.
    Strict mode raises an immediate exception.
    """
     
    # authors or moderators may edit a post
    if user.is_anonymous():
        writeable = False
    else:
        writeable = user.profile.can_moderate or (post.author == user)
    
    if strict and not writeable:
        msg = 'user %s edit access to post %s denied' % (user.id, post.id)
        logger.error(msg)
        raise Exception(msg)

    post.writeable = writeable
    return post.writeable

