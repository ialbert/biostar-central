from django import template
import urllib, hashlib


register = template.Library()

def comments(comments):
    return {'comments':comments}
    
def userlink(user):
    return {'user':user}
    
def userrep(user):
    return {'user':user}
    
def userbox(post, action='asked'):
    return {'post':post, 'action':action}
    
def gravatar(user, size=80):
    gravatar_url = "http://www.gravatar.com/avatar.php?"
    gravatar_url += urllib.urlencode({
        'gravatar_id':hashlib.md5(user.email).hexdigest(),
        'size':str(size)})
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, user.username)
    
register.inclusion_tag('widgets/comments.html')(comments)
register.inclusion_tag('widgets/userlink.html')(userlink)
register.inclusion_tag('widgets/userrep.html')(userrep)
register.inclusion_tag('widgets/userbox.html')(userbox)

register.simple_tag(gravatar)
