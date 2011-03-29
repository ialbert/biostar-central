from django import template
import urllib, hashlib
from biostar.server import models

register = template.Library()

def votebox(context, post):
    return { 'post':post,
            'upvoted': post.get_vote(context['user'], models.VOTE_UP) is not None,
            'downvoted':post.get_vote(context['user'], models.VOTE_DOWN) is not None}

def gravatar(user, size=80):
    gravatar_url = "http://www.gravatar.com/avatar.php?"
    gravatar_url += urllib.urlencode({
        'gravatar_id':hashlib.md5(user.email).hexdigest(),
        'size':str(size)})
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, user.username)
    
def comments(comments):
    return {'comments':comments}
    
def userlink(user):
    return {'user':user}
    
def userbox(post, action='asked'):
    return {'post':post, 'action':action}


register.simple_tag(gravatar)
register.inclusion_tag('widgets/votebox.html', takes_context=True)(votebox)
register.inclusion_tag('widgets/comments.html')(comments)
register.inclusion_tag('widgets/userlink.html')(userlink)
register.inclusion_tag('widgets/userbox.html')(userbox)




