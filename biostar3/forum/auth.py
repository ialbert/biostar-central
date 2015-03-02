"""
Access authorizations are performed here
"""
from string import strip
from .models import *
from django.utils.timezone import utc
from datetime import datetime

class AccessDenied(BaseException):
    pass

def now():
    return datetime.utcnow().replace(tzinfo=utc)

def tag_split(text):
    lower = lambda x: x.lower() if len(x) > 1 else x
    parts = text.split(",")
    parts = map(strip, parts)
    parts = filter(None, parts)
    parts = map(lower, parts)
    return parts

def create_toplevel_post(data, user, group):
    # Creating a top level post from  data
    title = data.get('title', '').strip()
    type = data.get('type', '')
    tags = data.get('tags', '')
    tags = tag_split(tags)
    content = data.get('content', '')

    # Create the post.
    post = Post.objects.create(content=content, title=title,
                              author=user, type=type, group=group)
    # Set the tags on the post
    post.tags.set(*tags)
    
    # Self referential ForeignKeys need to be updated explicitly!
    Post.objects.filter(pk=post.pk).update(root_id=post.id, parent_id=post.id)

    # Return the updated object, othewise the foreign keys are not set.
    post = Post.objects.get(pk=post.id)

    return post

def create_content_post(data, parent, post_type, user):
    # Creating a content level post from data
    content = data.get('content', '')
    post = Post.objects.create(parent=parent, content=content, type=post_type, author=user)
    return post

def read_access_post(user, post):
    """
    A user may read the post if the post is in a public group or
    the user is part of the group that the post was made in.
    """
    return post.root.group.public or user.usergroups.filter(name=post.root.group.name).exists()

def write_access_post(user, post):
    """
    A user may write the post if the post is readable and
    the user is the author of the post or is a moderator
    """
    write_cond = (user == post.author) or user.is_moderator
    return write_cond and read_access_post(user=user, post=post)

def thread_write_access(user, root):

    read_cond = read_access_post(post=root, user=user)

    def validator(user, post):
        write_cond = (user == post.author) or user.is_moderator
        return read_cond and write_cond

    return validator