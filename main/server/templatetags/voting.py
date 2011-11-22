from django import template
from main.server import const

register = template.Library()

def votebox(context, post, up, down):
    "This generates the up/down arrows"
    return { 'post':post, 
             'upvoted': post.id in up,
             'downvoted': post.id in down,
    }

register.inclusion_tag('widgets/vote.box.html', takes_context=True)(votebox)





