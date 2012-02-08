from django import template
from django.conf import settings
from django.template import Context, Template
from django.template.defaultfilters import stringfilter

register = template.Library()

@register.inclusion_tag('new/widgets/show.tags.html')
def show_tags(post):
    "Renders tags for a post"
    return {'post':post}

@register.inclusion_tag('new/widgets/post.info.small.html')
def post_info_small(action, post):
    "Renders post information"
    return {'post':post, 'action':action}
    
@register.inclusion_tag('new/widgets/tab.bar.html')
def tab_bar(request):
    "Renders post information"
    return { 'request': request }
    
@register.inclusion_tag('new/widgets/page.bar.html', takes_context=True)
def pagebar(context, anchor=''):
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