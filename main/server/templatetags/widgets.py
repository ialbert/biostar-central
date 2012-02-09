from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter

register = template.Library()

@register.inclusion_tag('new/widgets/show.tags.html')
def show_tags(post):
    "Renders tags for a post"
    return {'post':post}

@register.inclusion_tag('new/widgets/post.info.line.html')
def post_info_line(action, post):
    "Renders post information"
    return {'post':post, 'action':action}
    
@register.inclusion_tag('new/widgets/post.info.panel.html')
def post_info_panel(action, post):
    "Renders post information"
    return {'post':post, 'action':action}

@register.inclusion_tag('new/widgets/new.render.post.html', takes_context=True)
def render_post(context, post, tree):
    "Renders a post"
    return { 'post':post, 'tree':tree, 'request':context['request'] }

@register.inclusion_tag('new/widgets/render.child.html', takes_context=True)
def render_child(context, post, tree):
    "Renders a post"
    return { 'post':post, 'tree':tree, 'request':context['request'] }
   
    
@register.inclusion_tag('new/widgets/tab.bar.html')
def tab_bar(request):
    "Renders post information"
    return { 'request': request }
    
@register.inclusion_tag('new/widgets/page.bar.html', takes_context=True)
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