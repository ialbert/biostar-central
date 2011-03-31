from django import template
import urllib, hashlib
from datetime import datetime, timedelta


register = template.Library()

def comments(post):
    return {'post':post}
    
def userlink(user):
    return {'user':user}
    
def userrep(user):
    return {'user':user}
    
def userbox(post, action='asked'):
    return {'post':post, 'action':action}
    
def time_ago(time):
    delta = datetime.now() - time
    if delta < timedelta(minutes=1):
        return 'just now'
    if delta < timedelta(hours=1):
        return '%d min ago' % (delta.seconds // 60 )
    if delta < timedelta(days=1):
        return '%d hrs ago' % (delta.seconds // 3600 )
    if delta < timedelta(days=30):
        return '%d days ago' % delta.days
    return time.strftime('%b %d at %H:%M')
    
    
    
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
register.simple_tag(time_ago)

