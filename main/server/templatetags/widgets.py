from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter
from django.core.context_processors import csrf
import urllib
register = template.Library()


@register.simple_tag

@register.simple_tag
def show_value(value):
    return " (%s) " % value if value else ""
   
@register.simple_tag
def show_count(key, store):
    value = store.get(key)
    return "(%s)" % value if value else ""

@register.simple_tag
def show_type(post, flag):
    if flag:
        return "%s: " % post.get_type_display()
    else:
        return ""
   
@register.inclusion_tag('widgets/form.field.html',)
def form_field(field, label, help=''):
    errors = ", ".join(field.errors)
    return {'label':label, 'field':field, 'errors':errors, 'help':help}

@register.inclusion_tag('widgets/user.box.html',)
def user_box(user, size=100):
    return {'user':user, 'size':size}

@register.inclusion_tag('widgets/vote.box.html', takes_context=True)
def vote_box(context, post):
    return {'post':post}

@register.inclusion_tag('widgets/post.edit.actions.html')
def post_edit_actions(request, post):
    user = request.user
    user.can_moderate = user.is_authenticated() and user.profile.can_moderate
    return { 'post':post, 'request':request, 'user':user,  'SITE_DOMAIN':settings.SITE_DOMAIN}
    
@register.simple_tag
def search_url(term1, type='all'):
    patt  = '<a href="/search/?q=%s&t=%s">%s</a>'
    term2 = urllib.quote(term1)
    url = patt % (term2, type, term1)
    return url
    
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
    return { 'post':post, 'tree':tree, 'request':context['request'] }

@register.inclusion_tag('widgets/tab.bar.html')
def tab_bar(params={}, counts={}):
    "Renders the switchable tab on most pages"
    return { 'layout':params.get('layout'), 'tab': params.get('tab'), 'counts':counts, 'params':params }

@register.inclusion_tag('widgets/pill.bar.html')
def pill_bar(params={}, counts={}):
    "Renders the switchable pill bar on most pages"
    return { 'layout':params.get('layout'), 'pill':params.get('pill'), 'params':params, 'counts':counts }

@register.inclusion_tag('widgets/nav.bar.html', takes_context=True)
def nav_bar(context, user, params={}):
    "Renders top navigation bar"
    return { 'user':user, 'nav': params.get('nav'), 'q':params.get('q',''), 'request':context['request'] }

@register.inclusion_tag('widgets/page.bar.dropdown.html')
def page_bar_dropdown(selected, choices):
    "Renders top navigation bar"
    return { 'selected': selected, 'choices':choices }
    
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