from django import template
from django.conf import settings
import urllib, hashlib, re
from datetime import datetime, timedelta
from main.server import const, html, models, auth
from django.template import Context, Template, TemplateDoesNotExist
from django.core.context_processors import csrf
from urlparse import urlparse

register = template.Library()

from django.template.defaultfilters import stringfilter

uni = html.unicode_or_bust

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

@register.inclusion_tag('widgets/post.share.html')
def post_share(post):
    return {'post':post}

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
    
    username  = user.profile.display_name
    useremail = user.email.encode('utf8')

    gravatar_url = "http://www.gravatar.com/avatar.php?"
    gravatar_url += urllib.urlencode({
        'gravatar_id':hashlib.md5(useremail).hexdigest(),
        'size':str(size),
        'd':'identicon',
        }
    )
    return """<img src="%s" alt="gravatar for %s"/>""" % (gravatar_url, username)

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
        return 'Editor: '
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


templates = {}

def load_templates():
    for typeid, typename in const.POST_MAP.items():
        
        # this is the type of the template as a string
        typename = typename.lower()
        
        # see if the template has been overriden, and generate a default value
        default = 'rows/row.%s.html' % typename
        fname   = settings.TEMPLATE_ROWS.get(typename, default)
        try:
            templates[typeid] = template.loader.get_template(fname)
        except TemplateDoesNotExist:
            # fall back to a template that should exist
            #print "*** template loader loading default row for type '%s" % fname
            templates[typeid] = template.loader.get_template('rows/row.post.html')

# the template for the deleted row
row_deleted = template.loader.get_template('rows/row.deleted.html')

load_templates()

@register.simple_tag(takes_context=True)
def table_row(context, post, params, search_context=''):
    "Renders an html row for a post "
    global row_question, row_answer, row_comment, row_post, row_blog, row_forum, row_deleted

    row_deleted = template.loader.get_template('rows/row.deleted.html')

    if settings.DEBUG:
        # this is necessary to force the reload during development
        load_templates()

    c = Context( {"post": post, 'params':params, 'context': search_context, 'user':context['user']})
    if post.deleted:
        templ = row_deleted
    else:
        templ = templates[post.type]
    text = templ.render(c)
    return text