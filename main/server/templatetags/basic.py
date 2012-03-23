from django import template
from django.conf import settings
import urllib, hashlib, re
from datetime import datetime, timedelta
from main.server import const, html, models, auth
from django.template import Context, Template
from django.core.context_processors import csrf
from urlparse import urlparse

register = template.Library()

from django.template.defaultfilters import stringfilter

@register.filter(name='chunk')
@stringfilter
def quick_chunk(text, size=250):
    "Slices off text"
    return text[:size]

def smart_chunk(text):
    "Chunks by words"
    size, coll = 0, []
    for word in text.split():
        size += len(word)
        coll.append(word)
        if size > 180:
            break
    return ' '.join(coll)

@register.inclusion_tag('widgets/user.link.html')
def userlink(user):
    return {'user':user}
    
@register.inclusion_tag('widgets/tag.link.html')
def taglink(tag_name):
    return {'tag_name':tag_name}

@register.inclusion_tag('widgets/user.rep.html')
def userrep(user):
    return {'user':user}

@register.inclusion_tag('widgets/user.notes.html')
def usernotes(user):
    return { 'user':user }

@register.inclusion_tag('widgets/badge.icon.html')
def badgeicon(type):
    return {'type':type}

@register.inclusion_tag('widgets/action.box.html')
def actionbox(user, date, action='asked'):
    return {'user':user, 'date':date, 'action':action}
    
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
    if delta < timedelta(days=90):
        return '%d weeks ago' % int(delta.days/7)
    if delta < timedelta(days=730):
        return '%d months ago' % int(delta.days/30)
    # not quite exact
    diff = delta.days/365.0
    return '%0.1f years ago' % diff
    
    return time.strftime('%b %d at %H:%M')

@register.simple_tag
def gravatar(user, size=80):
    gravatar_url = "http://www.gravatar.com/avatar.php?"
    gravatar_url += urllib.urlencode({
        'gravatar_id':hashlib.md5(user.email).hexdigest(),
        'size':str(size),
        'd':'mm',
        }
    )
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, user.username)

@register.simple_tag(takes_context=True)
def navclass(context, ending):
    url  = context['request'].get_full_path()
    path = urlparse(url).path
    if  path == ending:
        return 'class="active"' 
    else:
        return ''

@register.simple_tag
def active(word, target):
    if word == target:
        return 'class="active"' 
    else:
        return ''
    
@register.simple_tag
def bignum(number):
    "Reformats numbers with qualifiers as K, M, G"
    try:
        value = float(number)/1000.0  
        if value > 10:
            return "%0.fk" % value
        elif value > 1:
            return "%0.1fk" % value
    except ValueError, exc:
        pass
    return str(number)
    
@register.simple_tag
def designation(user):
    "Renders a designation for the user"
    if user.profile.is_admin:
        return 'Administrator: '
    elif user.profile.is_moderator:
        return 'Moderator: '
    elif user.profile.type == const.USER_BLOG:
        return 'Blog: '
    return "User"

@register.simple_tag
def change_css(tag):
    "Changes django builtins to bootstrap classes"
    if tag == 'info' : return "alert-success"
    if tag == 'warning' : return "alert-info"
    if tag == 'error' : return "alert-error"
    return tag
    
@register.simple_tag
def flair(user):
    "Renders a designation for the user"
    if user.profile.is_admin:
        return '&diams;&diams;'
    elif user.profile.is_moderator:
        return '&diams;'
    return ""

# preload the templates 
row_post     = template.loader.get_template('rows/row.post.html')
row_answer   = template.loader.get_template('rows/row.answer.html')
row_comment  = template.loader.get_template('rows/row.comment.html')
row_blog     = template.loader.get_template('rows/row.blog.html')
row_question = template.loader.get_template('rows/row.question.html')

@register.simple_tag
def table_row(post):
    "Renders an html row for a post "
    global row_question, row_answer, row_comment
    
    if settings.DEBUG:
        # this is necessary to force the reload during development
        row_post     = template.loader.get_template('rows/row.post.html')
        row_answer   = template.loader.get_template('rows/row.answer.html')
        row_comment  = template.loader.get_template('rows/row.comment.html')
        row_blog     = template.loader.get_template('rows/row.blog.html')
        row_question = template.loader.get_template('rows/row.question.html')
   
    c = Context( {"post": post, 'root':post.root})
    if post.type == const.POST_BLOG:
        row = row_blog
    elif post.type == const.POST_QUESTION:
        row = row_question    
    elif post.type == const.POST_ANSWER:
        row = row_answer
    elif post.type == const.POST_COMMENT:
        row = row_comment
    else:
        row = row_post 
    
    text = row.render(c)
    return text