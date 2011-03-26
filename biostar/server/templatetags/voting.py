from django import template

register = template.Library()

def votebox(post):
    return { 'post':post }

register.inclusion_tag('widgets/votebox.html')(votebox)

