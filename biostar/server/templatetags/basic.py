from django import template
import urllib, hashlib
from datetime import datetime, timedelta

register = template.Library()

@register.inclusion_tag('widgets/comments.html')
def comments(post):
    return {'post':post}

@register.inclusion_tag('widgets/userlink.html')
def userlink(user):
    return {'user':user}
    
@register.inclusion_tag('widgets/taglink.html')
def taglink(tag_name):
    return {'tag_name':tag_name}

@register.inclusion_tag('widgets/userrep.html')
def userrep(user):
    return {'user':user}

@register.inclusion_tag('widgets/userbox.html')
def userbox(post, action='asked'):
    return {'post':post, 'action':action}

@register.simple_tag
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

@register.simple_tag
def gravatar(user, size=80):
    gravatar_url = "http://www.gravatar.com/avatar.php?"
    gravatar_url += urllib.urlencode({
        'gravatar_id':hashlib.md5(user.email).hexdigest(),
        'size':str(size)})
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, user.username)

@register.inclusion_tag('widgets/pagebar.html', takes_context=True)
def pagebar(context):
    return {
        'page': context['page'],
        'request': context['request'],
    }
    
@register.inclusion_tag('widgets/question-list-narrow.html')
def question_list_narrow(questions):
    return {'questions':questions}