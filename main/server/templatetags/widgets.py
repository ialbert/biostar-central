from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf

register = template.Library()

@register.simple_tag
def show_value(value):
    if value:
        return " (%s) " % value
    else:
        return ""
@register.inclusion_tag('widgets/form.field.html',)
def form_field(field, label, help=''):
    errors = ", ".join(field.errors)
    return {'label':label, 'field':field, 'errors':errors, 'help':help}

@register.inclusion_tag('widgets/user.box.html',)
def user_box(user, size=64):
    return {'user':user, 'size':size}

@register.inclusion_tag('widgets/vote.box.html', takes_context=True)
def vote_box(context, post):
    return {'post':post}

@register.inclusion_tag('widgets/post.edit.actions.html')
def post_edit_actions(request, post):
    return { 'post':post, 'request':request }
    
@register.inclusion_tag('widgets/show.tags.html')
def show_tags(post):
    "Renders tags for a post"
    return {'post':post}

@register.inclusion_tag('widgets/post.info.line.html')
def post_info_line(action, post):
    "Renders post information"
    return {'post':post, 'action':action}
    
@register.inclusion_tag('widgets/post.info.panel.html')
def post_info_panel(action, post):
    "Renders post information"
    return {'post':post, 'action':action}

@register.inclusion_tag('widgets/render.post.html', takes_context=True)
def render_post(context, post, tree):
    "Renders a post"
    return { 'post':post, 'tree':tree, 'request':context['request'], 'full':True }

@register.inclusion_tag('widgets/render.post.html', takes_context=True)
def render_answer(context, post, tree):
    "Renders a post"
    return { 'post':post, 'tree':tree, 'request':context['request'], 'full':False }

@register.inclusion_tag('widgets/tab.bar.html')
def tab_bar(tab=''):
    "Renders the switchable tab on most pages"
    return { 'tab': tab }

@register.inclusion_tag('widgets/nav.bar.html', takes_context=True)
def nav_bar(context, user, params={}):
    "Renders top navigation bar"
    return { 'user':user, 'nav': params.get('nav'), 'q':params.get('q',''), 'request':context['request'] }
    
@register.inclusion_tag('widgets/page.bar.html', takes_context=True)
def page_bar(context, anchor=''):
    path = context['request'].path
    return {
        'page'   : context['page'],
        'match'  : context.get('match',''),
        'query'  : context.get('query',''),
        'params' : context.get('params', ''),
        'request': context['request'],
        'anchor' : anchor,
        'path'   : path,
    }
    
# this contains the body of each comment
comment_body  = template.loader.get_template('widgets/render.comment.html')

@register.simple_tag
def comments(request, post, tree):
    global comment_body    
    if settings.DEBUG:
        # reload the template to get changes
        comment_body = template.loader.get_template('widgets/render.comment.html')
    if post.id in tree:
        text = render_comments(request=request, post=post, tree=tree)
    else:
        text = ''
    return text

def render_comments(request, post, tree):
    "Traverses the tree"
    global comment_body
            
    def traverse(node):
        out = [ '<div class="indent">' ]
        con = Context( {"post": node, 'user':request.user, 'request':request} )
        con.update(csrf(request))
        res = comment_body.render(con)
        out.append( res )
        for child in tree[node.id]:
            out.append( traverse(child) )        
        out.append( "</div>" )
        return '\n'.join(out)
    
    # this collects the comments for the post
    coll = []
    for node in tree[post.id]:
        coll.append( traverse(node) )
    return '\n'.join(coll)