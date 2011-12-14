from django import template
from django.conf import settings
import urllib, hashlib, re
from datetime import datetime, timedelta
from main.server import const, html, models, auth
from django.template import Context, Template
from django.core.context_processors import csrf

register = template.Library()

from django.template.defaultfilters import stringfilter

@register.filter(name='chunk')
@stringfilter
def quick_chunk(text):
    "Slices off text"
    return text[:250]

def smart_chunk(text):
    "Chunks by words"
    size, coll = 0, []
    for word in text.split():
        size += len(word)
        coll.append(word)
        if size > 180:
            break
    return ' '.join(coll)

@register.inclusion_tag('widgets/comments.html', takes_context=True)
def comments2(context, user, post):
    
    coll = []
    for comment in post.comments():
        comment.writeable = auth.authorize_post_edit(user=user, post=comment, strict=False)
        coll.append( comment )
       
    return { 'post':post, 'user':user, 'comments':coll}

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

@register.inclusion_tag('widgets/edit.box.html')
def editbox(request, post):
    return { 'post':post, 'request':request }
    
@register.inclusion_tag('widgets/badge.icon.html')
def badgeicon(type):
    return {'type':type}

@register.inclusion_tag('widgets/action.box.html')
def actionbox(user, date, action='asked'):
    return {'user':user, 'date':date, 'action':action}


  
@register.inclusion_tag('widgets/user.box.html')
def userbox(user):
    return {'user':user}
    
@register.inclusion_tag('widgets/revision.box.html')
def revisionbox(post):
    return {'post':post}

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

@register.inclusion_tag('bars/page.bar.html', takes_context=True)
def pagebar(context, anchor=''):
    path = context['request'].path
    return {
        'page': context['page'],
        'm': context.get('m',''),
        'q': context.get('q',''),
        'params': context.get('params', ''),
        'request': context['request'],
        'anchor':anchor,
        'path':path,
    }
    
@register.inclusion_tag('widgets/answer-list-narrow.html')
def answer_list_narrow(x):
    return {'answers':x}

@register.inclusion_tag('widgets/vote.box.html', takes_context=True)
def vote_box(context, post):
    return {'post':post}
    
@register.simple_tag(takes_context=True)
def navclass(context, include_path, exclude_paths=''):
    path = context['request'].get_full_path()
    if re.search(include_path, path):
        if not exclude_paths or (True not in [pat in path for pat in exclude_paths.split(' ')]):
            return 'class="youarehere"'
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
        return 'Administrator'
    elif user.profile.is_moderator:
        return 'Moderator'
    return "Registered user"
    
@register.simple_tag
def flair(user):
    "Renders a designation for the user"
    if user.profile.is_admin:
        return '&diams;&diams;'
    elif user.profile.is_moderator:
        return '&diams;'
    return ""

@register.inclusion_tag('widgets/render.post.html', takes_context=True)
def render_post(context, request, post, tree):
    return { 'post':post, 'tree':tree, 'request':request }
 
comment_body  = template.loader.get_template('widgets/comment.html')

def render_comments(request, post, tree):
    global comment_body
            
    def traverse(node):
        out = [ '<div class="indent">' ]
        con = Context( {"post": node, 'user':request.user} )
        con.update(csrf(request))
        res = comment_body.render(con)
        out.append( res )
        for child in tree[node.id]:
            out.append( traverse(child) )        
        out.append( "</div>" )
        return '\n'.join(out)
    
    # this collects only the comments
    coll = []
    for node in tree[post.id]:
        coll.append( traverse(node) )
    return '\n'.join(coll)

@register.simple_tag
def comments(request, post, tree):
    global comment_body    
    if settings.DEBUG:
        comment_body = template.loader.get_template('widgets/comment.html')
    if post.id in tree:
        text = render_comments(request=request, post=post, tree=tree)
    else:
        text = ''
    return text
    


# preload the templates 
row_question = template.loader.get_template('rows/row.question.html')
row_answer   = template.loader.get_template('rows/row.answer.html')
row_comment  = template.loader.get_template('rows/row.comment.html')

@register.simple_tag
def table_row(post):
    "Renders an html row for a post "
    global row_question, row_answer, row_comment
    
    if settings.DEBUG:
        # this is necessary to force the reload during development
        row_question = template.loader.get_template('rows/row.question.html')
        row_answer   = template.loader.get_template('rows/row.answer.html')
        row_comment  = template.loader.get_template('rows/row.comment.html')

    if post.type == const.POST_QUESTION:
        c = Context( {"post": post} )
        row = row_question
    elif post.type == const.POST_ANSWER:
        c = Context( {"post": post, 'root':post.root})
        row = row_answer
    else:
        c = Context( {"post": post, 'root':post.root})
        row = row_comment 
    
    text = row.render(c)
    return text